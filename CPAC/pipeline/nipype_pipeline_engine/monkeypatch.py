def patch_base_interface():
    """
    Monkey-patch nipype.interfaces.base.core.BaseInterface.run().

    This is a temporary workaround to prevent nipype from emitting
    .proc-* files in the current working directory.
    """
    from nipype.interfaces.base.core import (
        BaseInterface,
        config,
        indirectory,
        InterfaceResult,
        os,
        RuntimeContext,
        str2bool,
        write_provenance,
    )

    class PatchedBaseInterface(BaseInterface):
        def run(self, cwd=None, ignore_exception=None, **inputs):
            with indirectory(cwd or os.getcwd()):  # This is the only line that changed
                rtc = RuntimeContext(
                    resource_monitor=config.resource_monitor and self.resource_monitor,
                    ignore_exception=ignore_exception
                    if ignore_exception is not None
                    else self.ignore_exception,
                )

            with indirectory(cwd or os.getcwd()):
                self.inputs.trait_set(**inputs)
            self._check_mandatory_inputs()
            self._check_version_requirements(self.inputs)

            with rtc(self, cwd=cwd, redirect_x=self._redirect_x) as runtime:
                # Grab inputs now, as they should not change during execution
                inputs = self.inputs.get_traitsfree()
                outputs = None
                # Run interface
                runtime = self._pre_run_hook(runtime)
                runtime = self._run_interface(runtime)
                runtime = self._post_run_hook(runtime)
                # Collect outputs
                outputs = self.aggregate_outputs(runtime)

            results = InterfaceResult(
                self.__class__,
                rtc.runtime,
                inputs=inputs,
                outputs=outputs,
                provenance=None,
            )

            # Add provenance (if required)
            if str2bool(config.get("execution", "write_provenance", "false")):
                # Provenance will only throw a warning if something went wrong
                results.provenance = write_provenance(results)

            self._duecredit_cite()

            return results

    # Apply patch
    import nipype.interfaces.base.core as base_core
    base_core.BaseInterface.run = PatchedBaseInterface.run

