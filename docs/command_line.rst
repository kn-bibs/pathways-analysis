
**********************
Command Line Interface
**********************


User guide
==========

The command line interface has built in help. To display the help, please append `-h` to the program call, for example::

    ./patapy.py -h

The help option responds to arguments your provide, so you can get details about your method of choice with::

    ./patapy.py gsea -h

where gsea is the name of a method;
likewise, you can display help for any of samples specification options (case/control/data), e.g.::

    ./patapy.py control -h

Predefined parsers
==================

Parsers are defined in :mod:`command_line.main` module.

.. currentmodule:: command_line.main

.. automodule:: command_line.main
    :members:


Creating custom arguments and parsers
=====================================

Please use :mod:`command_line.parser` module to create custom parsers and arguments.

.. currentmodule:: command_line.parser

.. automodule:: command_line.parser
    :members:
