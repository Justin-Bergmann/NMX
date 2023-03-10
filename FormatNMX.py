from __future__ import annotations

from shutil import copy
from sys import argv
from typing import Tuple
import h5py
from scipy.constants import Planck, m_n

import cctbx.array_family.flex
from dials.array_family import flex


from dxtbx import IncorrectFormatError
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Detector, Goniometer, PolyBeam, Sequence






class FormatNMX(FormatHDF5):

    """
    Class to read files from NMX
    preprocessed files in scipp to abtain binned data
    """

    def __init__(self, image_file):
        if not FormatNMX.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        self.image_file = image_file
        self.nxs_file = self.open_file(image_file)
        self.detector = None
        self.raw_data = None      
        print("test3",self.nxs_file["NMX_data"].attrs["name"].encode('utf-8', 'surrogateescape').decode('latin-1'))


    
    def open_file(self, image_file):
        return h5py.File(image_file, "r")

    @staticmethod
    def understand(image_file):
        print("check NMX format1")
        try:
            print("understanding NMX format", FormatNMX.is_nmx_file(image_file) )
            return FormatNMX.is_nmx_file(image_file)
        except IOError:
            return False

    @staticmethod
    def is_nmx_file(image_file) -> bool:
        print("check NMX format")
        print("check2",image_file)

        """
        Confirms if image_file is from NMX
        """
        print("check NMX format")
        def get_name(image_file):
            print("checking NMX111")
     
            with h5py.File(image_file, "r") as handle:     
                print("handel1",handle["NMX_data"].attrs["name"].encode('utf-8', 'surrogateescape').decode('latin-1')) 
                     
                return handle["NMX_data"].attrs["name"].encode('utf-8', 'surrogateescape').decode('latin-1')
#            print("test read NXTOFRAW",FormatNXTOFRAW.understand(image_file))
#        if not FormatNXTOFRAW.understand(image_file):
#            return False
        return get_name(image_file) == "NMX"
    

    def get_spectra_idx_1D(self, panel: int, x_px: int, y_px: int) -> int:

        image_size = self._get_panel_size_in_px()
        total_pixels = image_size[0] * image_size[1]
        # Index offset in SXD data
        # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
        idx_offset = 361600
        panel_idx = (y_px * image_size[1]) + x_px
        panel_start_idx = (total_pixels * panel) + (idx_offset * panel)
        return int(panel_start_idx + panel_idx)



    def get_raw_spectra(self, normalize_by_proton_charge=True):
        """ loads TOF data for all detectores at once"""
        if normalize_by_proton_charge:
            proton_charge = self.nxs_file["NMX_data"]["proton_charge"][0]
            return (
                self.nxs_file["NMX_data"]["detector_1"]["counts"][:] / proton_charge
            )
        return self.nxs_file["NMX_data"]["detector_1"]["counts"][:]

    def save_spectra(self, spectra, output_filename):
        copy(self.image_file, output_filename)
        nxs_file = h5py.File(output_filename, "r+")
        del nxs_file["NMX_data"]["detector_1"]["counts"]
        nxs_file["NMX_data/detector_1"].create_dataset(
            "counts", spectra.shape, dtype=np.dtype("f8")
        )
        nxs_file["NMX_data/detector_1/counts"][:] = spectra
        nxs_file.close()    


    def load_raw_data(self, as_numpy_arrays=False, normalise_by_proton_charge=True):
        def get_detector_idx_array(detector_number, image_size, idx_offset):
            total_pixels = image_size[0] * image_size[1]
            min_range = (total_pixels * (detector_number - 1)) + (
                idx_offset * (detector_number - 1)
            )
            max_range = min_range + total_pixels
            return np.arange(min_range, max_range).reshape(image_size).T

        dataset = "NMX_data"
        raw_counts = self.nxs_file[dataset]["detector_1"]["counts"][0, :, :]

        if normalise_by_proton_charge:
            proton_charge = self.nxs_file[dataset]["proton_charge"][0]
            raw_counts = raw_counts / proton_charge

        num_panels = self._get_num_panels()
        image_size = self._get_panel_size_in_px()
        num_bins = len(self.get_tof_in_seconds())

        # Index offset in SXD data
        # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
        idx_offset = 4
        raw_data = []

        for n in range(1, num_panels + 1):
            idx_array = get_detector_idx_array(n, image_size, idx_offset)
            panel_array = np.zeros((idx_array.shape[0], idx_array.shape[1], num_bins))
            for c_i, i in enumerate(idx_array):
                for c_j, j in enumerate(i):
                    panel_array[c_i, c_j, :] = raw_counts[j, :]
            if as_numpy_arrays:
                raw_data.append(panel_array)
            else:
                flex_array = flex.double(np.ascontiguousarray(panel_array))
                flex_array.reshape(flex.grid(panel_array.shape))
                raw_data.append(flex_array)

        return tuple(raw_data)

    def get_raw_data(self, index):

        if self.raw_data is None:
            self.raw_data = self.load_raw_data()

        raw_data_idx = []
        for i in self.raw_data:
            arr = i[:, :, index : index + 1]
            arr.reshape(flex.grid(i.all()[0], i.all()[1]))
            raw_data_idx.append(arr)
        return tuple(raw_data_idx)

    def get_tof_in_seconds(self):
        return self._get_time_channels_in_seconds()

    def _get_time_channels_in_seconds(self):
        bins = self._get_time_channel_bins()
        return [(bins[i] + bins[i + 1]) * 0.5 * 10**-6 for i in range(len(bins) - 1)]

    def _get_time_channel_bins(self):
        return self.nxs_file["NMX_data"]["detector_1"]["tbin"][:]
    
    def _get_time_channels_in_seconds(self):
        bins = self._get_time_channel_bins()
        return [(bins[i] + bins[i + 1]) * 0.5 * 10**-6 for i in range(len(bins) - 1)]
    

    def _get_time_channels_in_usec(self):
        bins = self._get_time_channel_bins()
        return [(bins[i] + bins[i + 1]) * 0.5 for i in range(len(bins) - 1)]
    
    def get_tof_in_seconds(self):
        return self._get_time_channels_in_seconds()
    
    def get_tof_range(self):
        return (0, len(self.get_tof_in_seconds()))

    def get_wavelength_channels_in_ang(self):
        time_channels = self._get_time_channels_in_seconds()
        L = self._get_sample_to_moderator_distance() * 10**-3
        return [self.get_tof_wavelength_in_ang(L, i) for i in time_channels]
    
    def get_wavelength_channels(self):
        time_channels = self._get_time_channels_in_seconds()
        L = self._get_sample_to_moderator_distance() * 10**-3
        return [self.get_tof_wavelength(L, i) for i in time_channels]

    def get_wavelength_channels_in_A(self):
        wavelengths = self.get_wavelength_channels()
        return [i * 10**10 for i in wavelengths]
    
    def get_tof_wavelength(self, L, tof):
        return (Planck * tof) / (m_n * L)

    def get_tof_wavelength_in_ang(self, L, tof):
        return self.get_tof_wavelength(L, tof) * 10**10

    def get_max_slice_index(self):
        return len(self.get_tof_in_seconds()) - 1


    def get_max_slice_index(self):
        return len(self.get_tof_in_seconds()) - 1

    def get_pixel_spectra(self, panel_idx, x, y):
        if self.raw_data is None:
            self.raw_data = self.load_raw_data()

        time_channels = self._get_time_channels_in_usec()
        return time_channels, list(self.raw_data[panel_idx][x : x + 1, y : y + 1, :])

    def get_raw_spectra_two_theta(self):
        detector = self.get_detector()
        unit_s0 = self.get_beam().get_unit_s0()
        # Index offset in SXD data
        # See p24 of https://www.isis.stfc.ac.uk/Pages/sxd-user-guide6683.pdf
        idx_offset = 4
        raw_spectra_two_theta = []
        for panel in detector:
            two_theta = panel.get_two_theta_array(unit_s0)
            raw_spectra_two_theta += list(two_theta)
            raw_spectra_two_theta += [0] * idx_offset

        return raw_spectra_two_theta

    def get_detector(self, idx=None):
        return self._get_detector()
    
    def _get_detector(self):

        """
        Returns a  Detector instance with parameters taken from

        """

        num_panels = self._get_num_panels()
        panel_names = self._get_panel_names()
        panel_type = self._get_panel_type()
        image_size = self._get_panel_size_in_px()
        trusted_range = self._get_panel_trusted_range()
        pixel_size = self._get_panel_pixel_size_in_mm()
        fast_axes = self._get_panel_fast_axes()
        slow_axes = self._get_panel_slow_axes()
        panel_origins = self._get_panel_origins()
        gain = self._get_panel_gain()
        panel_projections = self._get_panel_projections_2d()
        detector = Detector()
        root = detector.hierarchy()

        for i in range(num_panels):
            panel = root.add_panel()
            panel.set_type(panel_type)
            panel.set_name(panel_names[i])
            panel.set_image_size(image_size)
            panel.set_trusted_range(trusted_range)
            panel.set_pixel_size(pixel_size)
            panel.set_local_frame(fast_axes[i], slow_axes[i], panel_origins[i])
            panel.set_gain(gain)
            r, t = panel_projections[i + 1]
            r = tuple(map(int, r))
            t = tuple(map(int, t))
            panel.set_projection_2d(r, t)

        return detector

    def _get_num_panels(self):
        #reads number of detector panales from file
            print("nr of detectors1: ",self.nxs_file['NMX_data/instrument'].attrs["nr_detector"], type(handle['NMX_data/instrument'].attrs["nr_detector"]))
            
            return self.nxs_file['NMX_data/instrument'].attrs["nr_detector"]

    def _get_panel_names(self):
        return ["%02d" % (i + 1) for i in range(self.nxs_file['NMX_data/instrument'].attrs["nr_detector"])]
    
    def _get_panel_type(self):
        return "SENSOR_PAD"
    


    def _get_panel_size_in_px(self):
        #number of pixels in each panal direction
        return (1280, 1280)

    def _get_panel_trusted_range(self):
        # 4 * 1280**2 plus buffer
        return (-1, 7000000)

    def _get_panel_pixel_size_in_mm(self):
        # pixes sice in mm
        return (0.4, 0.4)

    def _get_panel_fast_axes(self):
    #    has to be modified with variable orientation of detectors???
        return (
            self._panel_0_params()["fast_axis"],
            (-1.0, -0.0, -0.0),
            (-1.0, -0.0, -0.0),
            (-1.0, -0.0, -0.0),
        )

    def _get_panel_slow_axes(self):
    #    has to be modified with variable orientation of detectors  ???  
        return (
            self._panel_0_params()["slow_axis"],
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
   
        )

    def _get_panel_origins(self):
        #    has to be modified with variable orientation of detectors   ???
        return (
            self._panel_0_params()["origin"],
            (224.999, -96.0, 96.0),
            (60.809, -96.0, 236.945),
            (-214.172, -96.0, 118.198),
        )
    
    def _get_panel_projections_2d(self) -> dict:

        """
        Returns a projection of the
        detector flattened around the bottom panel (11), with
        all panels facing upright.

        """
#   adjustion to moving detecors
        p_w, p_h = self._get_panel_size_in_px()
        panel_pos = {
            0: ((1, 0, 0, 1), (0, 0)),
            1: ((1, 0, 0, 1), (p_h, 0)),
            2: ((1, 0, 0, 1), (-p_h, 0)),

        }

        return panel_pos

    def get_raw_data(self, index) -> Tuple:

        """
        Returns data of each panel for the ToF histogram bin at index
        """
        pass

    def get_beam(self, idx=None) -> PolyBeam:
        
        pass

    def get_beam(self, idx=None):
        sample_to_source_dir = self._get_sample_to_source_direction()
        sample_to_mod_d = self._get_sample_to_moderator_distance()
        wavelength_range = self._get_wavelength_range()
        return PolyBeamFactory.make_beam(
            sample_to_source_direction=sample_to_source_dir,
            sample_to_moderator_distance=sample_to_mod_d,
            wavelength_range=wavelength_range,
        )

    
    def _get_sample_to_source_direction(self):
        return (0, 0, -1)
    
    def _get_wavelength_range(self):
        '''wawelength range in A'''
        return (1.8, 3.55)
    
    def _get_sample_to_moderator_distance(self):
        """ gets distance between smaple and source in cm (moderator is not imlimentet jet ????)  """
        dist = abs(self.nxs_file['NMX_data/NXsource/distance'[:][0]])*100
        return dist

    def _get_panel_gain(self):
        return 1.0


    
    def _get_panel_size_in_mm(self):
        size_in_px = self._get_panel_size_in_px()
        pixel_size_in_mm = self._get_panel_pixel_size_in_mm()
        return tuple(
            [size_in_px[i] * pixel_size_in_mm[i] for i in range(len(size_in_px))]
        )

    def get_spectra_L1s(self, detector):

        pixel_size = self._get_panel_pixel_size_in_mm()
        panel_size = self._get_panel_size_in_px()
        num_pixels = panel_size[0] * panel_size[1]
        pixels = cctbx.array_family.flex.vec2_double(num_pixels)
        offset = (0, 0, 0, 0)
        count = 0
        for i in range(panel_size[0]):
            for j in range(panel_size[1]):
                pixels[count] = (i * pixel_size[0], j * pixel_size[1])

        spectra_L1 = []
        for p in range(self._get_num_panels()):
            s1 = detector[p].get_lab_coord(pixels)
            L1 = s1.norms() * 10**-3
            L1 = np.append(L1, offset)
            spectra_L1 = np.append(spectra_L1, L1)

        return spectra_L1

    def get_momentum_correction(self, detector):
        def get_momentum_bin_widths(L0, L1, tof_bins):
            wavelengths = [
                self.get_tof_wavelength_in_ang(L0 + L1, tof) for tof in tof_bins
            ]

            momentum = [2 * np.pi / i for i in wavelengths]
            bin_widths = np.abs(np.diff(momentum))
            return bin_widths

        spectra_L1 = self.get_spectra_L1s(detector)

        L0 = self._get_sample_to_moderator_distance() * 10**-3
        tof_bins = self._get_time_channel_bins()
        tof_bins = [i * 10**-6 for i in tof_bins]

        # For each pixel, get the bin width in 2pi/lambda
        correction = np.zeros((1, len(spectra_L1), len(tof_bins) - 1))
        for i, L1 in enumerate(spectra_L1):
            correction[0][i] = get_momentum_bin_widths(L0, L1, tof_bins)

        return correction

    def get_goniometer(self, idx=None) -> Goniometer:
        pass

#    def get_sequence(self, idx=None) -> Sequence:
#        pass
    def get_sequence(self, idx=None):
        image_range = (1, self.get_num_images())
        tof_in_seconds = self.get_tof_in_seconds()
        return SequenceFactory.make_tof_sequence(
            image_range=image_range,
            tof_in_seconds=tof_in_seconds,
            wavelengths=self.get_wavelength_channels_in_A(),
        )




if __name__ == "__main__":
    for arg in argv[1:]:
        print(FormatNMX.understand(arg))
