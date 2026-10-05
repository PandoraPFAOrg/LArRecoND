import numpy as np
import awkward as awk

import uproot as ur
import h5py as h5

from h5flow.data import dereference
import h5flow
import os
import sys
try:
    import yaml
except ImportError:
    yaml = None

# Refactored version of h5 to ROOT conversion, Bruce Howard - 2025
# initial scripts are now saved in e.g. _minirun6_2.py and _minirun6_3.py versions e.g. and thanks to Richie Diurba and any others who made these scripts
#
# Sep 2026 update: hit MERGING is now performed here (moved out of ndlar_flow by Sindhu). The
# converter reads the chosen *base* hits (prompt or filtered), merges same-channel
# hits within MERGE_CUT, and merges the per-hit truth backtracking (segment_ids +
# charge fraction) accordingly, so backtracking is preserved for merged hits.

# NOTE for more on uproot TTree writing see the below or the uproot documentation
# see https://stackoverflow.com/questions/72187937/writing-trees-number-of-baskets-and-compression-uproot

# ------------------------------------------------------------------ #
#  Merging configuration (defaults; overridden by the config yaml, arg 3)
# ------------------------------------------------------------------ #
# These are the Pandora-side defaults, used when no config yaml is
# passed. The merger was removed from ndlar_flow and handed to this step.
DEFAULT_BASE_HITS = 'filtered'   # 'prompt' or 'filtered'
DEFAULT_DO_MERGE  = True         # merge same-channel hits here
DEFAULT_MERGE_CUT = 65           # CRS ticks; MUST match the old ndlar_flow CalibHitMerger


def load_hit_config(path):
    '''Read the hit-selection/merging config yaml. Returns (base_hits, do_merge, merge_cut).
    Falls back to the module defaults for any missing key.'''
    cfg = {}
    if path and os.path.isfile(path):
        if yaml is None:
            raise RuntimeError('PyYAML not available but a config file was given: %s' % path)
        with open(path, 'r') as fh:
            cfg = yaml.safe_load(fh) or {}
    base = str(cfg.get('base_hits', DEFAULT_BASE_HITS)).lower()
    if base not in ('prompt', 'filtered'):
        base = DEFAULT_BASE_HITS
    return base, bool(cfg.get('merge', DEFAULT_DO_MERGE)), int(cfg.get('merge_cut', DEFAULT_MERGE_CUT))


def pick_base_key(f, prefer):
    '''Choose which base hit dataset to read.

    :param f: open h5py.File
    :param prefer: 'prompt' or 'filtered'
    :returns: one of 'prompt' / 'filtered' / 'final', falling back if the
              preferred dataset is absent. 'final' is kept for backwards
              compatibility with flow files produced before the
              calib_final_hits -> calib_filtered_hits rename.
    '''
    order = ['filtered', 'final', 'prompt'] if prefer == 'filtered' else ['prompt', 'filtered', 'final']
    for k in order:
        if 'charge/calib_%s_hits' % k in f:
            return k
    return 'prompt'


def merge_hits_and_backtrack(z, y, x, Q, E, ts, iog, ioc, chip, chan, ids,
                             bt_fraction=None, bt_segids=None, merge_cut=DEFAULT_MERGE_CUT):
    '''Merge same-channel hits (and their truth backtracking) within merge_cut.

    Hits are grouped by physical channel (io_group, z, y). Within a group,
    consecutive hits (sorted by ts) whose neighbour gap is < merge_cut are
    combined into one hit:
      - Q, E              -> summed
      - x, ts             -> charge-weighted mean
      - z, y, io_group,
        io_channel, chip_id, channel_id  -> channel-constant, taken as-is

    Backtracking (optional, MC only): each base hit carries a fixed-width list of
    (segment_ids, fraction) where fraction is that segment's share of the *base
    hit's* charge. For a merged hit the per-segment fraction is recomputed as the
    charge-weighted combination of its constituents' fractions, deduplicated by
    segment_id, kept up to the same number of slots (largest fraction first).

    Returns merged arrays in the same order plus (merged_fraction, merged_segids)
    (both None when no backtracking was supplied).
    '''
    n = len(z)
    if n == 0:
        return (z, y, x, Q, E, ts, iog, ioc, chip, chan, ids, bt_fraction, bt_segids)

    # ---- group base hits into merged hits ----
    order = np.lexsort((ts, y, z, iog))
    zc, yc, ic, tc = z[order], y[order], iog[order], ts[order].astype(np.int64)
    same = (zc[1:] == zc[:-1]) & (yc[1:] == yc[:-1]) & (ic[1:] == ic[:-1])
    close = np.abs(np.diff(tc)) < merge_cut
    starts = np.r_[True, ~(same & close)]         # a new merged hit begins here
    gsorted = np.cumsum(starts) - 1
    group = np.empty(n, dtype=np.int64)
    group[order] = gsorted
    ng = int(gsorted[-1]) + 1

    w = np.abs(Q).astype('float64')
    tot = np.zeros(ng, dtype='float64')
    np.add.at(tot, group, w)
    totw = np.where(tot == 0., 1e-300, tot)

    def _sum(a):
        s = np.zeros(ng, dtype='float64'); np.add.at(s, group, a.astype('float64')); return s

    def _wmean(a):
        s = np.zeros(ng, dtype='float64'); np.add.at(s, group, a.astype('float64') * w); return s / totw

    def _first(a):
        out = np.empty(ng, dtype=a.dtype)
        out[group[::-1]] = a[::-1]                 # reversed write -> lowest-index hit wins
        return out

    m_z = _first(z); m_y = _first(y)
    m_x = _wmean(x).astype('float32'); m_ts = _wmean(ts).astype('float32')
    m_Q = _sum(Q).astype('float32'); m_E = _sum(E).astype('float32')
    m_iog = _first(iog); m_ioc = _first(ioc); m_chip = _first(chip); m_chan = _first(chan)
    m_ids = np.arange(ng)

    # ---- merge the truth backtracking ----
    m_frac = m_seg = None
    if bt_fraction is not None and bt_segids is not None:
        n_slots = bt_fraction.shape[1]
        m_frac = np.zeros((ng, n_slots), dtype='float32')
        m_seg = np.full((ng, n_slots), -1, dtype=bt_segids.dtype)
        for g in range(ng):
            idx = np.where(group == g)[0]
            acc = {}                               # segment_id -> summed (fraction * charge)
            for i in idx:
                Qi = abs(float(Q[i]))
                fr = bt_fraction[i]; sd = bt_segids[i]
                for f_, s_ in zip(fr, sd):
                    if f_ != 0.:
                        acc[s_] = acc.get(s_, 0.) + float(f_) * Qi
            if not acc:
                continue
            denom = sum(abs(v) for v in acc.values()) or 1e-300
            top = sorted(acc.items(), key=lambda kv: -abs(kv[1]))[:n_slots]
            for j, (s_, c_) in enumerate(top):
                m_seg[g, j] = s_
                m_frac[g, j] = c_ / denom

    return (m_z, m_y, m_x, m_Q, m_E, m_ts, m_iog, m_ioc, m_chip, m_chan, m_ids, m_frac, m_seg)


# Main function with command line settable params
def printUsage():
    print('python h5_to_root_ndlarflow.py FileList IsData HitConfig LegacyMode OutName')
    print('-- Parameters')
    print('FileList    [REQUIRED]:                                         comma separated set of files to convert - note it will be one output')
    print('IsData      [OPTIONAL, DEFAULT = 0, is MC]:                     1 = Data, otherwise = MC')
    print('HitConfig   [OPTIONAL, DEFAULT = built-in filtered+merge]:      path to a hit-selection yaml (base_hits: prompt|filtered, merge: bool, merge_cut: int). Also accepts the strings "prompt"/"filtered" for convenience. See pandora_hits_config.yaml.')
    print('LegacyMode  [OPTIONAL, DEFAULT = 0, no legacy]:                 0 = no legacy mode, 1 = samples < MiniRun6, 2 = > MiniRun6 but no usec time')
    print('OutName     [OPTIONAL, DEFAULT = input[0]+"_hits_uproot.root"]: string for an output file name if you want to override. Note that default writes to current directory.')
    print('')
    print('NOTE: The output of this file should then be processed with the rootToRootConversion macro to get the format expected by LArRecoND.')
    print('')

def main(argv=None):
    fileNames=[]
    useData=False
    basePref=DEFAULT_BASE_HITS   # preferred base hits: 'prompt' or 'filtered'
    do_merge=DEFAULT_DO_MERGE
    merge_cut=DEFAULT_MERGE_CUT
    legacyMode=0
    overrideOutname=1
    outname=''

    MeV2GeV=0.001
    trueXOffset=0 # Offsets if geometry changes
    trueYOffset=0 #42+268
    trueZOffset=0 #-1300

    if len(sys.argv)==1:
        print('---------------------------------------------------------------')
        print('Must at least pass a file location/name to be converted, usage:')
        print('---------------------------------------------------------------')
        printUsage()
        return
    if len(sys.argv)>1:
        if str(sys.argv[1])=='help' or str(sys.argv[1])=='h' or str(sys.argv[1])=='-h' or str(sys.argv[1])=='--help':
            print('---------------------------------------------------------------')
            print('usage:')
            print('---------------------------------------------------------------')
            printUsage()
            return
        elif sys.argv[1]!=None:
            fileList=str(sys.argv[1])
            fileNames=fileList.split(',')
        if len(sys.argv)>2 and sys.argv[2]!=None:
            if int(sys.argv[2])==1:
                useData=True
        if len(sys.argv)>3 and sys.argv[3]!=None:
            arg3 = str(sys.argv[3])
            if arg3 in ('prompt', 'filtered'):
                basePref = arg3
            elif os.path.isfile(arg3):
                basePref, do_merge, merge_cut = load_hit_config(arg3)
            else:
                # back-compat: 0 = prompt base, otherwise = filtered base
                try:
                    basePref = 'prompt' if int(arg3)==0 else 'filtered'
                except ValueError:
                    print('Could not interpret HitConfig arg "%s"; using defaults (%s, merge=%s, cut=%d)'
                          % (arg3, basePref, do_merge, merge_cut))
        if len(sys.argv)>4 and sys.argv[4]!=None:
            legacyMode=int(sys.argv[4])
        if len(sys.argv)>5 and sys.argv[5]!=None:
            outname=str(sys.argv[5])
            overrideOutname=0

    MaxArrayDepth=int(10000)
    MaxArrayDepthData=int(100000)
    isWritten=False

    if overrideOutname==1:
        outname = fileNames[0].split('/')[-1]+'_hits_uproot.root'

    ## We are choosing to write a bogus subevent to set the types of all the branches.
    ## The hope is this will then work well even in the case where the first event we'd see is actually a bad event
    ##################################################
    if len(fileNames) > 0:
        # Simple versions of the input vectors where everything is set to 0 of the proper type
        hits_z = np.array([0.]).astype('float32')
        hits_y = np.array([0.]).astype('float32')
        hits_x = np.array([0.]).astype('float32')
        hits_Q = np.array([0.]).astype('float32')
        hits_E = np.array([0.]).astype('float32')
        hits_ts = np.array([0.]).astype('float32')
        hits_io_group = np.array([0.]).astype('uint8')
        hits_io_channel = np.array([0.]).astype('uint8')
        hits_chip_id = np.array([0.]).astype('uint8')
        hits_channel_id = np.array([0.]).astype('uint8')
        runID = np.array( [0], dtype='int32' )
        subrunID = np.array( [0], dtype='int32' )
        eventID = np.array( [0], dtype='int32' )
        triggerID = np.array( [0], dtype='int32')
        event_start_t = np.array( [-5], dtype='int32' )
        event_end_t = np.array( [-5], dtype='int32' )
        event_unix_ts = np.array( [-5], dtype='int32' )
        if legacyMode!=1 and legacyMode!=2:
            event_unix_ts_usec = np.array( [-5], dtype='int32' )
        if useData==False:
            matches = np.array( [0] ).astype('uint16')
            packetFrac = np.array( [0.] ).astype('float32')
            pdgHit = np.array( [0] ).astype('int32')
            trackID = np.array( [0] ).astype('int64')
            particleID = np.array( [0] ).astype('int64')
            particleIDLocal = np.array( [0] ).astype('int64')
            interactionIndex = np.array( [0] ).astype('int64')
            trajStartX = np.array( [0.] ).astype('float32')
            trajStartY = np.array( [0.] ).astype('float32')
            trajStartZ = np.array( [0.] ).astype('float32')
            trajEndX = np.array( [0.] ).astype('float32')
            trajEndY = np.array( [0.] ).astype('float32')
            trajEndZ = np.array( [0.] ).astype('float32')
            trajID = np.array( [0] ).astype('int64')
            trajIDLocal = np.array( [0] ).astype('int64')
            trajPDG = np.array( [0] ).astype('int32')
            trajE = np.array( [0.] ).astype('float32')
            trajPx = np.array( [0.] ).astype('float32')
            trajPy = np.array( [0.] ).astype('float32')
            trajPz = np.array( [0.] ).astype('float32')
            trajVertexID = np.array( [0] ).astype('int64')
            trajParentID = np.array( [0] ).astype('int64')
            nu_vtx_id = np.array([0]).astype('int64')
            nu_vtx_x = np.array([0.]).astype('float32')
            nu_vtx_y = np.array([0.]).astype('float32')
            nu_vtx_z = np.array([0.]).astype('float32')
            nu_vtx_E = np.array([0.]).astype('float32')
            nu_pdg = np.array([0]).astype('int32')
            nu_px = np.array([0.]).astype('float32')
            nu_py = np.array([0.]).astype('float32')
            nu_pz = np.array([0.]).astype('float32')
            nu_iscc = np.array([0]).astype('int32')
            nu_code = np.array([0]).astype('int32')

        # Set up the dictionaries to write to the file
        event_dict = { 'run':runID, 'subrun':subrunID, 'event':eventID, "triggers":triggerID, 'unix_ts':event_unix_ts,
                       'event_start_t':event_start_t, 'event_end_t':event_end_t }
        if legacyMode!=1 and legacyMode!=2:
            event_dict['unix_ts_usec'] = event_unix_ts_usec

        if useData==False:
            other_dict = {  'x':hits_x, 'y':hits_y, 'z':hits_z, 'ts':hits_ts, 'io_group':hits_io_group, 'io_channel':hits_io_channel , 'chip_id':hits_chip_id, 'channel_id':hits_channel_id, 'charge':hits_Q, 'E':hits_E, 'matches':matches,\
                            'mcp_energy':trajE, 'mcp_pdg':trajPDG, 'mcp_nuid':trajVertexID, 'mcp_vertex_id':trajVertexID,\
                            'mcp_idLocal':trajIDLocal, 'mcp_id':trajID, 'mcp_px':trajPx, 'mcp_py':trajPy, 'mcp_pz':trajPz,\
                            'mcp_mother':trajParentID, 'mcp_startx':trajStartX, 'mcp_starty':trajStartY, 'mcp_startz':trajStartZ,\
                            'mcp_endx':trajEndX, 'mcp_endy':trajEndY, 'mcp_endz':trajEndZ,\
                            'nuID':nu_vtx_id, 'vertex_id':nu_vtx_id, 'nue':nu_vtx_E, 'nuPDG':nu_pdg,\
                            'nupx':nu_px, 'nupy':nu_py, 'nupz':nu_pz, 'nuvtxx':nu_vtx_x, 'nuvtxy':nu_vtx_y,\
                            'nuvtxz':nu_vtx_z, 'mode':nu_code, 'ccnc':nu_iscc,\
                            'hit_packetFrac':packetFrac, 'hit_particleID':particleID, 'hit_particleIDLocal':particleIDLocal,\
                            'hit_pdg':pdgHit, 'hit_vertexID':interactionIndex, 'hit_segmentID':trackID }
        else:
            other_dict = {  'x':hits_x, 'y':hits_y, 'z':hits_z, 'ts':hits_ts, 'io_group':hits_io_group, 'io_channel':hits_io_channel, 'chip_id':hits_chip_id, 'channel_id':hits_channel_id, 'charge':hits_Q, 'E':hits_E }

        max_entries=0
        for key in other_dict.keys():
            if len(other_dict[key]) > max_entries:
                max_entries = len(other_dict[key])

        if useData==True:
            nSubEvents = int(max_entries/MaxArrayDepthData)+1
            for idxSubEvent in range(nSubEvents):
                first = MaxArrayDepth*idxSubEvent
                last = MaxArrayDepth*(idxSubEvent+1)
                event_dict['subevent'] = np.array([idxSubEvent], dtype='int32')
                for key in other_dict.keys():
                    event_dict[key] = awk.values_astype(awk.Array([other_dict[key][first:last]]),other_dict[key].dtype)
                fout = ur.recreate(outname)
                fout['subevents'] = event_dict
                isWritten=True
        else:
            nSubEvents = int(max_entries/MaxArrayDepth)+1
            for idxSubEvent in range(nSubEvents):
                first = MaxArrayDepth*idxSubEvent
                last = MaxArrayDepth*(idxSubEvent+1)
                event_dict['subevent'] = np.array([idxSubEvent], dtype='int32')
                for key in other_dict.keys():
                    event_dict[key] = awk.values_astype(awk.Array([other_dict[key][first:last]]),other_dict[key].dtype)
                fout = ur.recreate(outname)
                fout['subevents'] = event_dict
                isWritten=True
            del packetFrac
            del particleID
            del particleIDLocal
            del pdgHit
            del interactionIndex
            del trackID
    ##################################################

    for fileIdx in range(len(fileNames)):
        print('Processing file',fileIdx,'of',len(fileNames))
        fileName=fileNames[fileIdx]

        f = h5.File(fileName)
        events=f['charge/events/data']
        flow_out=h5flow.data.H5FlowDataManager(fileName,"r")

        # Decide which base hits to read for this file (prompt/filtered, with fallback)
        baseKey = pick_base_key(f, basePref)
        print('Using base hits: charge/calib_%s_hits  (merge=%s, merge_cut=%d)' % (baseKey, do_merge, merge_cut))

        eventsToRun=len(events)

        # Get the array of the trigger type for every event in the file
        triggerIDsData=flow_out["charge/events","charge/ext_trigs",events["id"][:]]
        triggerIDsAll=np.array(np.ma.getdata(triggerIDsData["iogroup"]),dtype='int32')

        triggerIDs = np.array( np.broadcast_to( (1 << 31) - 1, shape=eventsToRun ) )
        for i in range(len(triggerIDsAll)):
            if np.sum(triggerIDsAll[i])==0:
                continue
            if 5 in triggerIDsAll[i]:
                triggerIDs[i] = 5
            else:
                triggerIDs[i] = triggerIDsAll[i][0]

        for ievt in range(eventsToRun):
            badEvt=False

            if ievt%10==0:
                print('Currently on',ievt,'of',eventsToRun)
            event = events[ievt]
            event_base_hits=flow_out["charge/events/","charge/calib_"+baseKey+"_hits", events["id"][ievt]]

            if len(event_base_hits[0])==0:
                print('This event seems empty in the hits array, setting as bad event. Trigger type (',triggerIDs[ievt],')')
                badEvt=True

            # Read the base hits for this event
            #######################################
            if badEvt==False:
                # Check if the only values are masked and call this a bad event if so
                if np.ma.count_masked(event_base_hits["z"][0]) == len(event_base_hits[0]):
                    print('This event has a hit z array ( len hits =', len(event_base_hits[0]), ') that appears to be only masked values, setting as bad event. Trigger type (',triggerIDs[ievt],')')
                    badEvt=True

            if badEvt==False:
                hits_z = (np.ma.getdata(event_base_hits["z"][0])+trueZOffset).astype('float32')
                hits_y = ( np.ma.getdata(event_base_hits["y"][0])+trueYOffset ).astype('float32')
                hits_x = ( np.ma.getdata(event_base_hits["x"][0])+trueXOffset ).astype('float32')
                hits_Q = ( np.ma.getdata(event_base_hits["Q"][0]) ).astype('float32')
                hits_E = ( np.ma.getdata(event_base_hits["E"][0]) ).astype('float32')
                hits_ts = ( np.ma.getdata(event_base_hits["ts_pps"][0]) ).astype('float32')
                hits_io_group = ( np.ma.getdata(event_base_hits["io_group"][0]) ).astype('uint8')
                hits_io_channel = ( np.ma.getdata(event_base_hits["io_channel"][0]) ).astype('uint8')
                hits_chip_id = ( np.ma.getdata(event_base_hits["chip_id"][0]) ).astype('uint8')
                hits_channel_id = ( np.ma.getdata(event_base_hits["channel_id"][0]) ).astype('uint8')
                hits_ids = np.ma.getdata(event_base_hits["id"][0])
            else:
                hits_z = np.array([]).astype('float32')
                hits_y = np.array([]).astype('float32')
                hits_x = np.array([]).astype('float32')
                hits_Q = np.array([]).astype('float32')
                hits_E = np.array([]).astype('float32')
                hits_ts = np.array([]).astype('float32')
                hits_io_group = np.array([]).astype('uint8')
                hits_io_channel = np.array([]).astype('uint8')
                hits_chip_id = np.array([]).astype('uint8')
                hits_channel_id = np.array([]).astype('uint8')
                hits_ids = np.array([])

            if badEvt==False and len(hits_ids)<2:
                print('This event has < 2 hit IDs, setting as bad event. Trigger type (',triggerIDs[ievt],')')
                badEvt=True

            # Note: in a few places below, we will want to have handy the spillID
            spillID = 0
            if useData==False and badEvt==False:
                unmaskedSpillIDs = []
                if baseKey=='prompt':
                    allSpillIDs=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hits_ids]["event_id"]
                    unmaskedSpillIDs = allSpillIDs.data[ ~allSpillIDs.mask ]
                else:
                    event_hits_prompt=flow_out["charge/events/","charge/calib_prompt_hits", events["id"][ievt]]
                    hits_ids_prompt = np.ma.getdata(event_hits_prompt["id"][0])
                    allSpillIDs=flow_out["charge/calib_prompt_hits","charge/packets","mc_truth/segments",hits_ids_prompt]["event_id"]
                    unmaskedSpillIDs = allSpillIDs.data[ ~allSpillIDs.mask ]
                if len(unmaskedSpillIDs) > 0:
                    spillID = unmaskedSpillIDs[0]
                else:
                    print('This event has no spillID from matches that we want to use in grabbing true particles/neutrinos. Setting as bad event. Trigger type (',triggerIDs[ievt],')')
                    badEvt=True

            # ---- read base backtracking (MC), BEFORE merging (uses base hit ids) ----
            base_frac = base_seg = None
            if useData==False and badEvt==False:
                base_charge_path = 'charge/calib_%s_hits' % baseKey
                base_truth_path  = 'mc_truth/calib_%s_hit_backtrack' % baseKey
                base_bt = flow_out[base_charge_path, base_truth_path, hits_ids[:]][:,0]
                base_frac = np.ma.getdata(base_bt['fraction'])
                base_seg  = np.ma.getdata(base_bt['segment_ids'])

            # ---- merge hits (and backtracking) ----
            merged_frac = merged_seg = None
            if badEvt==False and do_merge:
                (hits_z, hits_y, hits_x, hits_Q, hits_E, hits_ts,
                 hits_io_group, hits_io_channel, hits_chip_id, hits_channel_id, hits_ids,
                 merged_frac, merged_seg) = merge_hits_and_backtrack(
                    hits_z, hits_y, hits_x, hits_Q, hits_E, hits_ts,
                    hits_io_group, hits_io_channel, hits_chip_id, hits_channel_id, hits_ids,
                    base_frac, base_seg, merge_cut)
            else:
                merged_frac, merged_seg = base_frac, base_seg

            # Start with the non-spill info, this is all ~like the current form
            #   but not repeating
            runID = np.array( [0], dtype='int32' )
            subrunID = np.array( [0], dtype='int32' )
            eventID = np.array( [event['id']], dtype='int32' )
            triggerID = np.array( [triggerIDs[ievt]], dtype='int32')

            maxTimeFromTrigger = np.max( np.array( [event['ts_start']], dtype='float64' ) )
            useTimeFromTrigger = True
            if maxTimeFromTrigger > ((1 << 31)-1):
                useTimeFromTrigger = False

            if (triggerID!=((1 << 31)-1)) | (useTimeFromTrigger==True):
                event_start_t = np.array( [event['ts_start']], dtype='int32' )
                event_end_t = np.array( [event['ts_end']], dtype='int32' )
                event_unix_ts = np.array( [event['unix_ts']], dtype='int32' )
                if legacyMode!=1 and legacyMode!=2:
                    event_unix_ts_usec = np.array( [event['unix_ts_usec']], dtype='int32' )
            else:
                event_start_t = np.array( [-5], dtype='int32' )
                event_end_t = np.array( [-5], dtype='int32' )
                event_unix_ts = np.array( [-5], dtype='int32' )
                if legacyMode!=1 and legacyMode!=2:
                    event_unix_ts_usec = np.array( [-5], dtype='int32' )

            if useData==False:
                # Truth-level info for hits (from the MERGED backtracking)
                #######################################
                if badEvt==False and merged_frac is not None:
                    # Matches (per merged hit) + flattened contributing-segment lists
                    nonzero = (merged_frac != 0.)
                    matches = nonzero.sum(axis=1).astype('uint16')
                    packetFrac = merged_frac[nonzero].astype('float32')
                    segmentIDs = merged_seg[nonzero]

                    all_segments = f['mc_truth/segments/data']
                    all_segments = all_segments[ np.where(all_segments['event_id']==spillID) ]
                    all_segmentIDs = all_segments['segment_id']

                    lookup = {}
                    for i, v in enumerate(all_segmentIDs):
                        if v not in lookup:
                            lookup[v] = i

                    segments_where = np.array([lookup[s] for s in segmentIDs if s in lookup], dtype='int64')

                    pdgHit = all_segments['pdg_id'][segments_where].astype('int32')
                    trackID = all_segments['segment_id'][segments_where].astype('int64')
                    particleID = all_segments['file_traj_id'][segments_where].astype('int64')
                    particleIDLocal = all_segments['traj_id'][segments_where].astype('int64')
                    interactionIndex = all_segments['vertex_id'][segments_where].astype('int64')
                else:
                    matches = np.zeros( len(hits_z), dtype='uint16' )
                    packetFrac = np.array( [] ).astype('float32')
                    pdgHit = np.array( [] ).astype('int32')
                    trackID = np.array( [] ).astype('int64')
                    particleID = np.array( [] ).astype('int64')
                    particleIDLocal = np.array( [] ).astype('int64')
                    interactionIndex = np.array( [] ).astype('int64')

                # Truth-level info for the spill
                #######################################
                if badEvt==False:
                    # Trajectories
                    traj_indicesArray = np.where(flow_out['mc_truth/trajectories/data']["event_id"] == spillID)[0]
                    traj = flow_out["mc_truth/trajectories/data"][traj_indicesArray]
                    trajStartX = (traj['xyz_start'][:,0]).astype('float32')
                    trajStartY = (traj['xyz_start'][:,1]).astype('float32')
                    trajStartZ = (traj['xyz_start'][:,2]).astype('float32')
                    trajEndX = (traj['xyz_end'][:,0]).astype('float32')
                    trajEndY = (traj['xyz_end'][:,1]).astype('float32')
                    trajEndZ = (traj['xyz_end'][:,2]).astype('float32')
                    trajID = (traj['file_traj_id']).astype('int64')
                    trajIDLocal = (traj['traj_id']).astype('int64')
                    trajPDG = (traj['pdg_id']).astype('int32')
                    trajE = (traj['E_start']*MeV2GeV).astype('float32')
                    trajPx = (traj['pxyz_start'][:,0]*MeV2GeV).astype('float32')
                    trajPy = (traj['pxyz_start'][:,1]*MeV2GeV).astype('float32')
                    trajPz = (traj['pxyz_start'][:,2]*MeV2GeV).astype('float32')
                    trajVertexID = (traj['vertex_id']).astype('int64')
                    trajParentID = (traj['parent_id']).astype('int64')
                else:
                    trajStartX = np.array( [] ).astype('float32')
                    trajStartY = np.array( [] ).astype('float32')
                    trajStartZ = np.array( [] ).astype('float32')
                    trajEndX = np.array( [] ).astype('float32')
                    trajEndY = np.array( [] ).astype('float32')
                    trajEndZ = np.array( [] ).astype('float32')
                    trajID = np.array( [] ).astype('int64')
                    trajIDLocal = np.array( [] ).astype('int64')
                    trajPDG = np.array( [] ).astype('int32')
                    trajE = np.array( [] ).astype('float32')
                    trajPx = np.array( [] ).astype('float32')
                    trajPy = np.array( [] ).astype('float32')
                    trajPz = np.array( [] ).astype('float32')
                    trajVertexID = np.array( [] ).astype('int64')
                    trajParentID = np.array( [] ).astype('int64')

                # Vertices
                if badEvt==False:
                    vertex_indicesArray = np.where(flow_out["/mc_truth/interactions/data"]["event_id"] == spillID)[0]
                    vtx = flow_out["/mc_truth/interactions/data"][vertex_indicesArray]
                    nu_vtx_id = (vtx['vertex_id']).astype('int64')
                    if legacyMode!=1:
                        nu_vtx_x = (vtx['x_vert']).astype('float32')
                        nu_vtx_y = (vtx['y_vert']).astype('float32')
                        nu_vtx_z = (vtx['z_vert']).astype('float32')
                    else:
                        nu_vtx_x = (vtx['vertex'][:,0]).astype('float32')
                        nu_vtx_y = (vtx['vertex'][:,1]).astype('float32')
                        nu_vtx_z = (vtx['vertex'][:,2]).astype('float32')
                    nu_vtx_E = (vtx['Enu']*MeV2GeV).astype('float32')
                    nu_pdg = (vtx['nu_pdg']).astype('int32')
                    nu_px = (vtx['nu_4mom'][:,0]*MeV2GeV).astype('float32')
                    nu_py = (vtx['nu_4mom'][:,1]*MeV2GeV).astype('float32')
                    nu_pz = (vtx['nu_4mom'][:,2]*MeV2GeV).astype('float32')
                    # Little bit of gymnastics here
                    ccnc = vtx['isCC']
                    nu_iscc = np.invert(ccnc).astype('int32')
                    # And more gymnastics here
                    codes = 1000*np.ones(len(nu_vtx_id),dtype='int32')
                    idxQE = np.where(vtx['isQES']==True)
                    idxRES = np.where(vtx['isRES']==True)
                    idxDIS = np.where(vtx['isDIS']==True)
                    idxMEC = np.where(vtx['isMEC']==True)
                    idxCOH = np.where(vtx['isCOH']==True)
                    idxCOHQE = np.where((vtx['isCOH']==True) & (vtx['isQES']==True))
                    codes[idxQE] = 0
                    codes[idxRES] = 1
                    codes[idxDIS] = 2
                    codes[idxCOH] = 3
                    codes[idxCOHQE] = 4
                    codes[idxMEC] = 10
                    nu_code = codes
                else:
                    nu_vtx_id = np.array([]).astype('int64')
                    nu_vtx_x = np.array([]).astype('float32')
                    nu_vtx_y = np.array([]).astype('float32')
                    nu_vtx_z = np.array([]).astype('float32')
                    nu_vtx_E = np.array([]).astype('float32')
                    nu_pdg = np.array([]).astype('int32')
                    nu_px = np.array([]).astype('float32')
                    nu_py = np.array([]).astype('float32')
                    nu_pz = np.array([]).astype('float32')
                    nu_iscc = np.array([]).astype('int32')
                    nu_code = np.array([]).astype('int32')

            ## Rebuild now with all the individual types
            event_dict = { 'run':runID, 'subrun':subrunID, 'event':eventID, "triggers":triggerID, 'unix_ts':event_unix_ts,
                           'event_start_t':event_start_t, 'event_end_t':event_end_t }
            if legacyMode!=1 and legacyMode!=2:
                event_dict['unix_ts_usec'] = event_unix_ts_usec

            if useData==False:
                other_dict = {  'x':hits_x, 'y':hits_y, 'z':hits_z, 'ts':hits_ts, 'io_group':hits_io_group, 'io_channel':hits_io_channel,'chip_id':hits_chip_id , 'channel_id':hits_channel_id ,'charge':hits_Q, 'E':hits_E, 'matches':matches,\
                                'mcp_energy':trajE, 'mcp_pdg':trajPDG, 'mcp_nuid':trajVertexID, 'mcp_vertex_id':trajVertexID,\
                                'mcp_idLocal':trajIDLocal, 'mcp_id':trajID, 'mcp_px':trajPx, 'mcp_py':trajPy, 'mcp_pz':trajPz,\
                                'mcp_mother':trajParentID, 'mcp_startx':trajStartX, 'mcp_starty':trajStartY, 'mcp_startz':trajStartZ,\
                                'mcp_endx':trajEndX, 'mcp_endy':trajEndY, 'mcp_endz':trajEndZ,\
                                'nuID':nu_vtx_id, 'vertex_id':nu_vtx_id, 'nue':nu_vtx_E, 'nuPDG':nu_pdg,\
                                'nupx':nu_px, 'nupy':nu_py, 'nupz':nu_pz, 'nuvtxx':nu_vtx_x, 'nuvtxy':nu_vtx_y,\
                                'nuvtxz':nu_vtx_z, 'mode':nu_code, 'ccnc':nu_iscc,\
                                'hit_packetFrac':packetFrac, 'hit_particleID':particleID, 'hit_particleIDLocal':particleIDLocal,\
                                'hit_pdg':pdgHit, 'hit_vertexID':interactionIndex, 'hit_segmentID':trackID }
            else:
                other_dict = {  'x':hits_x, 'y':hits_y, 'z':hits_z, 'ts':hits_ts, 'io_group':hits_io_group, 'io_channel':hits_io_channel, 'chip_id':hits_chip_id, 'channel_id':hits_channel_id, 'charge':hits_Q, 'E':hits_E }

            max_entries=0
            for key in other_dict.keys():
                if len(other_dict[key]) > max_entries:
                    max_entries = len(other_dict[key])

            if useData==True:
                nSubEvents = int(max_entries/MaxArrayDepthData)+1
                for idxSubEvent in range(nSubEvents):
                    first = MaxArrayDepth*idxSubEvent
                    last = MaxArrayDepth*(idxSubEvent+1)
                    event_dict['subevent'] = np.array([idxSubEvent], dtype='int32')
                    for key in other_dict.keys():
                        event_dict[key] = awk.values_astype(awk.Array([other_dict[key][first:last]]),other_dict[key].dtype)
                    if isWritten==False:
                        print('TAKE NOTE! I thought I should have already made the output file by now, but I have "isWritten" as False, so I am attempting to create the output file.')
                        fout = ur.recreate(outname)
                        fout['subevents'] = event_dict
                        isWritten=True
                    else:
                        fout['subevents'].extend(event_dict)
            else:
                nSubEvents = int(max_entries/MaxArrayDepth)+1
                for idxSubEvent in range(nSubEvents):
                    first = MaxArrayDepth*idxSubEvent
                    last = MaxArrayDepth*(idxSubEvent+1)
                    event_dict['subevent'] = np.array([idxSubEvent], dtype='int32')
                    for key in other_dict.keys():
                        event_dict[key] = awk.values_astype(awk.Array([other_dict[key][first:last]]),other_dict[key].dtype)
                    if isWritten==False:
                        print('TAKE NOTE! I thought I should have already made the output file by now, but I have "isWritten" as False, so I am attempting to create the output file.')
                        fout = ur.recreate(outname)
                        fout['subevents'] = event_dict
                        isWritten=True
                    else:
                        fout['subevents'].extend(event_dict)
                del packetFrac
                del particleID
                del particleIDLocal
                del pdgHit
                del interactionIndex
                del trackID

        fout.close()
        print('end of code')

if __name__=="__main__":
    main()
