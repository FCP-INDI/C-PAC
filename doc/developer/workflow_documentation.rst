.. _workflow_documentation:

**********************
Workflow Documentation
**********************
The bulk of documentation for workflows is to be done in the source file itself.  Such practice will minimize the number of documentation files and ensure documentation is kept upto date with code changes.


.. _documenting_workflows:

How to write workflow code with proper documentation
====================================================

Proper workflow documentation will follow the general standards provided by numpy (https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).

:doc:`Example of a documented workflow <../workflows/anatpreproc>`

**Required Information:**

*Parameters*: refer to workflow inputs

*Returns*: refer to workflow outputs.

*Notes*: should provide an in-depth description of the proccessing which occurs inside a workflow.

*Examples*: should provide at least one example of how to appropriately use the workflow.  The example should be completely standalone and executable directly from a python interpreter with CPAC correctly installed.  The outputs of a workflow should be displayed to the user in a relevant manner.  For instance, the names of the files created or the numerical values of the resulting output run for example data.

Generating website documenation for workflows
=============================================

After workflow code has been properly documented, the documentation must be integrated into the browser viewable site.  Nearly all workflow documentation is generated directly from the source file.  This is made possible by the `autodoc extension of sphinx <http://sphinx.pocoo.org/ext/autodoc.html>`_ which converts the docstrings located in the source to html.

To add a documented workflow to the CPAC documentation site:

* In ./CPAC/doc/workflows/, create a *workflow_name*.rst file with the following text (`example <https://raw.github.com/openconnectome/C-PAC/master/doc/workflows/anatpreproc.rst>`_)::

    Human readable name of workflow (eg. Anatomical Preprocessing)
    ==============================================================
    
    .. currentmodule:: workflow_name
    .. autofunction:: create_workflow_name


* Add your workflow to the table of contents in ./CPAC/doc/workflows/index.rst (`example <https://raw.github.com/openconnectome/C-PAC/master/doc/workflows/index.rst>`_)::

    Workflows
    =========

    Contents:

    .. toctree::
       :maxdepth: 2
       ...
       workflow_name

* Make your workflow code available to sphinx by adding the following line to ./CPAC/doc/conf.py (`example <https://github.com/openconnectome/C-PAC/blob/master/doc/conf.py#L29>`_)::

    sys.path.append(os.path.abspath('sphinxext'))
    ...
    sys.path.append(os.path.abspath('../CPAC/workflow_name/'))

* Generate the html documentation using sphinx by running the command in ./CPAC/doc/::

    make html
