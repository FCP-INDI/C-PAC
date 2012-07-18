.. _workflow_documentation:

**********************
Workflow Documentation
**********************
The bulk of documentation for workflows is to be done in the source file itself.  Such practice will minimize the number of documentation files and ensure documentation is kept upto date with code changes.


.. _documenting_workflows:

How to write workflow code with proper documentation
====================================================

Proper workflow documentation will follow the general standards provided by numpy (https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).  Because workflow creation functions are themselves python functions, their parameters and what they return should be documented as such.  Please refer to the `example.py <https://github.com/numpy/numpy/blob/master/doc/example.py>`_ provided by numpy for syntax.


:doc:`Example of a documented workflow <../workflows/anat_preproc>`


**Notes section:**

Workflows have nipype inputs and outputs which are specific to the workflow returned.  Any workflow specific information should be documented in the *Notes* section.  At a minimum, workflows are required to have the following information:

* A link pointing to the file on GitHub.
* *Workflow Inputs*: refer to workflow inputs. This is often received by an IdentityInterface node named 'inputspec'.
* *Workflow Outputs*: refer to workflow outputs.  This is often received by an IdentityInterface node named 'outputspec'.
* Explanation of the processing which occurs inside the workflow

Example of *Notes* section for a workflow::

    Notes
    -----
    
    `Source <GitHub link pointing to workflow source file>`_
    
    Workflow Inputs::
    
        inputspec.X : possible data type(s) of X
            Detailed description of X
    
    Workflow Outputs::
    
        outputspec.Y : possible data type(s) of Y
            Detailed description of Y

    X is denoised and spatially smoothed.

**Examples Section:**

This section should provide at least one example of how to appropriately use the workflow.  The example should be completely standalone and executable directly from a python interpreter with CPAC correctly installed.  The outputs of a workflow should be displayed to the user in a relevant manner.  For instance, the names of the files created or the numerical values of the resulting output run for example data.

Generating website documenation for workflows
=============================================

After workflow code has been properly documented, the documentation must be integrated into the browser viewable site.  Nearly all workflow documentation is generated directly from the source file.  This is made possible by the `autodoc extension of sphinx <http://sphinx.pocoo.org/ext/autodoc.html>`_ which converts the docstrings located in the source to html.

To add a documented workflow to the CPAC documentation site:

* In ./CPAC/doc/workflows/, create a *workflow_name*.rst file with the following text (`example <https://raw.github.com/openconnectome/C-PAC/master/doc/workflows/anat_preproc.rst>`_)::

    Human readable name of workflow (eg. Anatomical Preprocessing)
    ==============================================================
    
    .. automodule:: CPAC.workflow_name
        :members:

* Add your workflow to the table of contents in ./CPAC/doc/workflows/index.rst (`example <https://raw.github.com/openconnectome/C-PAC/master/doc/workflows/index.rst>`_)::

    Workflows
    =========

    Contents:

    .. toctree::
       :maxdepth: 2
       ...
       workflow_name

* Generate the html documentation using sphinx by running the command in ./CPAC/doc/::

    make html
