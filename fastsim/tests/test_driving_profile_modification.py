"""
Tests that check the drive cycle modification functionality.
"""
import unittest

import numpy as np

import fastsim


MPH_TO_MPS = 1.0 / fastsim.params.mphPerMps


class TestDrivingProfileModification(unittest.TestCase):
    def setUp(self) -> None:
        # create a trapezoidal trip shape
        trapz = fastsim.cycle.make_cycle(
            [0.0, 10.0, 45.0, 55.0, 100.0],
            [0.0, 40.0 * MPH_TO_MPS, 40.0 * MPH_TO_MPS, 0.0, 0.0])
        trapz = fastsim.cycle.resample(trapz, new_dt=1.0)
        self.trapz = fastsim.cycle.Cycle(cyc_dict=trapz)
        self.veh = fastsim.vehicle.Vehicle(5)
        self.sim_drive = fastsim.simdrive.SimDriveClassic(self.trapz, self.veh)
        self.sim_drive_coast = fastsim.simdrive.SimDriveClassic(self.trapz, self.veh)
        return super().setUp()
    
    def tearDown(self) -> None:
        return super().tearDown()
    
    def test_that_impose_coast_causes_coast(self):
        "Test the standard interface to Eco-Approach for 'free coasting'"
        coast_start_mph = 39.99
        self.assertFalse(self.sim_drive.impose_coast.any(), "All impose_coast starts out False")
        while self.sim_drive.i < len(self.trapz.cycSecs):
            i = self.sim_drive.i
            prev_i = max(0, i-1)
            self.sim_drive_coast.impose_coast[i] = self.sim_drive_coast.impose_coast[prev_i] or self.sim_drive_coast.mphAch[prev_i] >= coast_start_mph
            self.sim_drive.sim_drive_step()
            self.sim_drive_coast.sim_drive_step()
        max_trace_miss_m__s = np.absolute(self.trapz.cycMps - self.sim_drive.mpsAch).max()
        max_trace_miss_coast_m__s = np.absolute(self.trapz.cycMps - self.sim_drive_coast.mpsAch).max()
        self.assertTrue(max_trace_miss_m__s < 0.01, f"Max trace miss: {max_trace_miss_m__s} m/s")
        self.assertTrue(max_trace_miss_coast_m__s > 1.0, f"Max trace miss: {max_trace_miss_coast_m__s} m/s")
        self.assertFalse(self.sim_drive.impose_coast.any())
        self.assertFalse(self.sim_drive_coast.impose_coast[0])
        self.assertTrue(self.sim_drive_coast.impose_coast[-1])
        if False:
            import matplotlib.pyplot as plt
            (fig, ax) = plt.subplots()
            ax.plot(self.trapz.cycSecs, self.trapz.cycMph, 'gray', label='shadow-trace')
            ax.plot(self.trapz.cycSecs, self.sim_drive_coast.mphAch, 'blue', label='coast')
            ax.plot(self.trapz.cycSecs, self.sim_drive_coast.mphAch, 'r.')
            ax.set_xlabel('Elapsed Time (s)')
            ax.set_ylabel('Speed (mph)')
            ax.legend()
            fig.tight_layout()
            fig.savefig('coasting-by-time.png', dpi=300)
            (fig, ax) = plt.subplots()
            ax.plot(self.trapz.cycDistMeters.cumsum(), self.trapz.cycMph, 'gray', label='shadow-trace')
            ax.plot(self.sim_drive_coast.distMeters.cumsum(), self.sim_drive_coast.mphAch, 'blue', label='coast')
            ax.plot(self.sim_drive_coast.distMeters.cumsum(), self.sim_drive_coast.mphAch, 'r.')
            ax.set_title('Coasting over Trapezoidal Cycle')
            ax.set_xlabel('Distance Traveled (m)')
            ax.set_ylabel('Speed (mph)')
            ax.legend()
            fig.tight_layout()
            fig.savefig('coasting-by-distance.png', dpi=300)
            fig = None
            ax = None
            (fig, ax) = plt.subplots()
            ax.plot(self.trapz.cycSecs, self.trapz.cycDistMeters.cumsum(), 'gray', label='shadow-trace')
            ax.plot(self.sim_drive_coast.cyc.cycSecs, self.sim_drive_coast.distMeters.cumsum(), 'blue', label='coast')
            ax.plot(self.sim_drive_coast.cyc.cycSecs, self.sim_drive_coast.distMeters.cumsum(), 'r.')
            ax.set_title('Coasting over Trapezoidal Cycle')
            ax.set_xlabel('Elapsed Time (s)')
            ax.set_ylabel('Distance Traveled (m)')
            ax.legend()
            fig.tight_layout()
            fig.savefig('coasting-distance-by-time.png', dpi=300)

    def test_eco_approach_modeling(self):
        "Test a simplified model of eco-approach"
        fastsim.simdrive.run_eco_approach(self.sim_drive_coast)
        self.assertFalse(self.sim_drive_coast.impose_coast.all(), "Assert we are not always in coast")
        self.assertTrue(self.sim_drive_coast.impose_coast.any(), "Assert we are at least sometimes in coast")
        max_trace_miss_coast_m__s = np.absolute(self.trapz.cycMps - self.sim_drive_coast.mpsAch).max()
        self.assertTrue(max_trace_miss_coast_m__s > 1.0, "Assert we deviate from the shadow trace")
        self.assertTrue(self.sim_drive_coast.mphAch.max() > 20.0, "Assert we at least reach 20 mph")
        #self.assertAlmostEqual(
        #    self.trapz.cycDistMeters.sum(),
        #    self.sim_drive_coast.distMeters.sum(), 1,
        #    "Assert the end distances are equal\n" +
        #    f"Got {self.trapz.cycDistMeters.sum()} m and {self.sim_drive_coast.distMeters.sum()} m")

