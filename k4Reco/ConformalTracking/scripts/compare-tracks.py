# A simple script to compare the tracks from the Gaudi and Marlin output
import argparse
from podio.root_io import Reader

parser = argparse.ArgumentParser(description="Compare tracks from Gaudi and Marlin")
parser.add_argument(
    "--gaudi-file", default="output_conformal_tracking.root", help="Gaudi output file"
)
parser.add_argument(
    "--marlin-file", default="output_REC.edm4hep.root", help="Marlin output file"
)

parser.add_argument(
    "--gaudi-tracks", default="NewSiTracks", help="Gaudi tracks collection"
)
parser.add_argument(
    "--marlin-tracks", default="SiTracks", help="Marlin tracks collection"
)

args = parser.parse_args()

reader_gaudi = Reader(args.gaudi_file)
reader_marlin = Reader(args.marlin_file)

events_gaudi = reader_gaudi.get("events")
events_marlin = reader_marlin.get("events")

for i, frame_gaudi in enumerate(events_gaudi):
    frame_marlin = events_marlin[i]
    tracks_gaudi = frame_gaudi.get(args.gaudi_tracks)
    tracks_marlin = frame_marlin.get(args.marlin_tracks)
    print(f"Checking event {i} with {len(tracks_gaudi)} tracks")
    assert len(tracks_gaudi) == len(
        tracks_marlin
    ), f"Number of tracks differ: {len(tracks_gaudi)} vs {len(tracks_marlin)}"
    for j, (track_gaudi, track_marlin) in enumerate(zip(tracks_gaudi, tracks_marlin)):
        print(f"Checking track {j}")
        assert (
            track_gaudi.getType() == track_marlin.getType()
        ), f"Type differ for track {j}: {track_gaudi.getType()} vs {track_marlin.getType()}"
        assert (
            track_gaudi.getChi2() == track_marlin.getChi2()
        ), f"Chi2 differ for track {j}: {track_gaudi.getChi2()} vs {track_marlin.getChi2()}"
        assert (
            track_gaudi.getNdf() == track_marlin.getNdf()
        ), f"Ndf differ for track {j}: {track_gaudi.getNdf()} vs {track_marlin.getNdf()}"
        assert (
            track_gaudi.getNholes() == track_marlin.getNholes()
        ), f"Nholes differ for track {j}: {track_gaudi.getNholes()} vs {track_marlin.getNholes()}"
        assert len(track_gaudi.getTrackerHits()) == len(
            track_marlin.getTrackerHits()
        ), f"Number of hits differ for track {j}: {len(track_gaudi.getTrackerHits())} vs {len(track_marlin.getTrackerHits())}"
        assert len(track_gaudi.getTrackStates()) == len(
            track_marlin.getTrackStates()
        ), f"Number of states differ for track {j}: {len(track_gaudi.getTrackStates())} vs {len(track_marlin.getTrackStates())}"
        for z, (hit_gaudi, hit_marlin) in enumerate(
            zip(track_gaudi.getTrackerHits(), track_marlin.getTrackerHits())
        ):
            assert (
                hit_gaudi.id() == hit_marlin.id()
            ), f"Hit {z} differ for track {j}: {hit_gaudi.id()} vs {hit_marlin.id()}"
        for z, (ts_gaudi, ts_marlin) in enumerate(
            zip(track_gaudi.getTrackStates(), track_marlin.getTrackStates())
        ):
            assert (
                ts_gaudi.location == ts_marlin.location
            ), f"Location differ for track {j}: {ts_gaudi.location} vs {ts_marlin.location}"
            assert (
                ts_gaudi.D0 == ts_marlin.D0
            ), f"D0 differ for track {j}: {ts_gaudi.D0} vs {ts_marlin.D0}"
            assert (
                ts_gaudi.phi == ts_marlin.phi
            ), f"Phi differ for track {j}: {ts_gaudi.phi} vs {ts_marlin.phi}"
            assert (
                ts_gaudi.omega == ts_marlin.omega
            ), f"Omega differ for track {j}: {ts_gaudi.omega} vs {ts_marlin.omega}"
            assert (
                ts_gaudi.Z0 == ts_marlin.Z0
            ), f"Z0 differ for track {j}: {ts_gaudi.Z0} vs {ts_marlin.Z0}"
            assert (
                ts_gaudi.tanLambda == ts_marlin.tanLambda
            ), f"TanLambda differ for track {j}: {ts_gaudi.tanLambda} vs {ts_marlin.tanLambda}"
            assert (
                ts_gaudi.referencePoint == ts_marlin.referencePoint
            ), f"ReferencePoint differ for track {j}: {ts_gaudi.referencePoint} vs {ts_marlin.referencePoint}"
            assert (
                ts_gaudi.covMatrix == ts_marlin.covMatrix
            ), f"CovMatrix differ for track {j}: {ts_gaudi.covMatrix} vs {ts_marlin.covMatrix}"
            # Time is set to -1 in the converter
            # assert ts_gaudi.time == ts_marlin.time, f"Time differ for track {j}: {ts_gaudi.time} vs {ts_marlin.time}"
