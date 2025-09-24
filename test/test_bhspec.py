import os
import json
import unittest
import urllib.request
import pycbc
from pycbc import conversions
from pycbc import inference
from pycbc.inference import io
from pycbc.inference.models import read_from_config
from pycbc.workflow import WorkflowConfigParser

class TestBHSpecModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # set the working directory
        whereami = os.path.dirname(os.path.realpath(__file__))
        cls.working_dir = f"{whereami}/../examples/inference/time_marg_bhspec"
        # the config file to use
        cls.config_file = "gw150914_hierarchical_gated.ini"
        # the expected values file
        cls.expected_values_file = "expected_maxl-gw150914.json"
        # url to the frame files
        cls.frame_files_url = \
            "https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/"
        # swith to the working directory
        os.chdir(cls.working_dir)
        # Load the configuration
        cls.cp = WorkflowConfigParser([cls.config_file])

        # Download frame files specified in the [ringdown__data] section
        frame_files = cls.cp.get("ringdown__data", "frame-files").split()
        for frame_file in frame_files:
            # strip off the instrument name
            frame_file = frame_file.split(":")[-1]
            if not os.path.exists(frame_file):
                url = os.path.join(cls.frame_files_url, frame_file)
                urllib.request.urlretrieve(url, frame_file)

        # Load expected parameter values and expected loglikelihood from
        # the JSON file
        with open(cls.expected_values_file, "r") as f:
            cls.expected_values = json.load(f)

        # turn on logging
        pycbc.init_logging(True)
        # Load the model
        # we'll load it with a copy because the original may be modified by
        # the model
        cls.model = read_from_config(cls.cp.__deepcopy__(cls.cp))


    def test_loglikelihood(self):
        """Tests that the loglikelihood and pre/post merger SNRs are close to
        the expected values."""
        # Set model parameters
        params = {p: self.expected_values[p]
                  for p in self.model.variable_params}
        self.model.update(**params)

        # Compute loglikelihood
        computed_loglikelihood = self.model.loglikelihood

        # Check if the computed loglikelihood is close to the expected value
        expected_loglikelihood = self.expected_values["loglikelihood"]
        self.assertAlmostEqual(
            computed_loglikelihood, expected_loglikelihood, places=2,
            msg=f"Loglikelihood mismatch: expected {expected_loglikelihood}, "
                f"got {computed_loglikelihood}"
        )

        # check the pre- and post-merger snrs
        postloglr = self.model.submodels['ringdown'].loglr
        preloglr = self.model.submodels['inspiral'].loglr

        # calculate inspiral and ringdown snrs
        snr_post = conversions.snr_from_loglr(postloglr)
        snr_pre = conversions.snr_from_loglr(preloglr)
        # check if the snrs are close to the expected values
        expected_snr_post = self.expected_values["ringdown__snr"]
        expected_snr_pre = self.expected_values["inspiral__snr"]
        self.assertAlmostEqual(
            snr_post, expected_snr_post, places=2,
            msg=f"Post-merger SNR mismatch: expected {expected_snr_post}, "
                f"got {snr_post}"
        )
        self.assertAlmostEqual(
            snr_pre, expected_snr_pre, places=2,
            msg=f"Pre-merger SNR mismatch: expected {expected_snr_pre}, "
                f"got {snr_pre}"
        )

    def test_normalization(self):
        """Tests that the model is properly normalized."""
        params = {p: self.expected_values[p]
                  for p in self.model.variable_params}
        self.model.update(**params)
        logl = self.model.loglikelihood
        # check that the loglikelihood is properly normalized if normalize
        # is set to True
        lognorm = 0
        for submodel in self.model.submodels.values():
            submodel.normalize = True
            for det in submodel.detectors:
                start_index, end_index = submodel.gate_indices(det)
                lognorm += submodel.det_lognorm(det, start_index, end_index)
        # check that the normalization is not zero
        self.assertNotEqual(lognorm, 0,
                            msg="Normalization is zero, check the model.")
        # now call the normed loglikelihood and check that it is the same
        # as the unnormed loglikelihood + lognorm
        self.model.update(**params)
        normed_logl = self.model.loglikelihood
        self.assertAlmostEqual(
            normed_logl, logl + lognorm,
            places=2,
            msg=f"Loglikelihood mismatch with normalization: expected "
                f"{logl + lognorm}, got {normed_logl}")

    def test_file_io(self):
        """Tests that we can write out a checkpoint file, then can load a
        model using the config, data, and psds from that checkpoint file and
        get the same loglikelihood."""
        # Write out a checkpoint file
        output_file = "test_output.hdf"
        checkpoint_file = f'{output_file}.checkpoint'
        if os.path.exists(checkpoint_file):
            os.remove(checkpoint_file)

        # Create a checkpoint file the way pycbc_inference would
        sampler = inference.sampler.load_from_config(
            self.cp, self.model, output_file=output_file, nprocesses=1)

        # make sure the checkpoint file has the expected name
        self.assertEqual(sampler.checkpoint_file, checkpoint_file)

        # add the config file
        with io.loadfile(checkpoint_file, 'a') as fp:
            fp.write_config_file(self.cp)

        # Now load back from the checkpoint file
        with io.loadfile(checkpoint_file, "r") as fp:
            # Load the data, psds, and config file from the checkpoint
            inspdata = fp.read_data(group='inspiral')
            rddata = fp.read_data(group='ringdown')
            insppsds = fp.read_psds(group='inspiral')
            rdpsds = fp.read_psds(group='ringdown')
            cp = fp.read_config_file()

        # Initialize a new model using the loaded data, psds, and config
        new_model = read_from_config(cp, inspiral__data=inspdata,
                                     inspiral__psds=insppsds,
                                     ringdown__data=rddata,
                                     ringdown__psds=rdpsds)

        # also check if the normalization is set in the config file, the loaded
        # model has normalization set
        cp.set("inspiral__model", "normalize", "")
        cp.set("ringdown__model", "normalize", "")
        normed_model = read_from_config(cp, inspiral__data=inspdata,
                                        inspiral__psds=insppsds,
                                        ringdown__data=rddata,
                                        ringdown__psds=rdpsds)
        # check that the normalization is True in the normed model
        for submodel in normed_model.submodels.values():
            self.assertTrue(submodel.normalize,
                            msg=f"Normalization not set in {submodel.name}")
        # Set parameters for the new models
        params = {p: self.expected_values[p]
                  for p in new_model.variable_params}
        new_model.update(**params)
        normed_model.update(**params)

        # Compute loglikelihood for the new model
        new_loglikelihood = new_model.loglikelihood
        normed_logl = normed_model.loglikelihood

        # Compare the loglikelihood of the new model with the original model
        self.model.update(**params)
        original_loglikelihood = self.model.loglikelihood
        self.assertAlmostEqual(
            new_loglikelihood, original_loglikelihood, places=2,
            msg=f"Loglikelihood mismatch after loading from checkpoint: "
                f"expected {original_loglikelihood}, got {new_loglikelihood}"
        )

        # Calculate the normalization of the normed model and compare to the
        # original model
        lognorm = 0
        for submodel in normed_model.submodels.values():
            submodel.normalize = True
            for det in submodel.detectors:
                start_index, end_index = submodel.gate_indices(det)
                lognorm += submodel.det_lognorm(det, start_index, end_index)
        # check that the normalization is not zero
        self.assertNotEqual(lognorm, 0,
                            msg="Normalization is zero, check the model.")
        # check that the loglikelihood is the same as the unnormed loglikelihood
        self.assertAlmostEqual(
            normed_logl - lognorm, new_loglikelihood,
            places=2,
            msg=f"Loglikelihood mismatch with normalization: expected "
                f"{new_loglikelihood}, got {normed_logl - lognorm}")

        # Clean up the checkpoint file
        for fn in [sampler.checkpoint_file, sampler.backup_file]:
            if os.path.exists(fn):
                os.remove(fn)

if __name__ == "__main__":
    unittest.main()