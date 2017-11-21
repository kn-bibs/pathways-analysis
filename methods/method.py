from abc import abstractmethod, ABC
from models import Experiment
from utils import AbstractRegisteringType, abstract_property


class MethodResult(ABC):
    """Result should contain list of matched pathways or processes

    to be displayed in the results table. The result can include
    additional information for the end user (in `description` field).

    The names of properties of the items (pathways or processes) in
    the list which should be used for table creation ought to be
    enlisted in `columns` property.

    Additional results created by a method should be presented in
    `files` field of the result.
    """

    @abstract_property
    def columns(self):
        pass

    def __init__(self, scored_list, files=None, description=''):
        self.scored_list = scored_list
        self.files = files or []
        self.description = description


class Method(metaclass=AbstractRegisteringType):
    """Defines method of pathway analysis & its arguments.

    Simple arguments (like ``threshold``) can be simply defined as
    arguments and keyword arguments of `__init__`.

    For example::

        class MyMethod(Method)
            def __init__(threshold:float=0.05):
                pass

    For the simple arguments following information will be deduced:
        - type: will be retrieved from type annotations; currently only
          non-abstract types (int, str, float and so on) are supported.
          We can implement abstract types from `typing` if needed.
        - default: from keyword arguments.
        - help: will be retrieved from docstrings

    If you need more advanced options (like aggregation), or just do not
    like having a mess in your `__init__` signature, please define the
    arguments in body of your class using :class:`~command_line.parser.Argument` constructor.

    For example::

        class MyMethod(Method):

            database = Argument(
                type=argparse.FileType('r'),
                help='Path to file with the database'
            )

            def __init__(threshold:float=0.05, database=None):
                pass

    If help is given in both :class:`~command_line.parser.Argument` and docstring,
    then the help from `Argument()` takes precedence over the help in docstrings
    (as docstrings should cover not only CLI usage but also describe how to use
    the method as a standalone object - to enable advanced users to customize methods).
    """

    @abstract_property
    def help(self):
        """Return string providing help for this method.

        The help message shows up when `./run method_name -h`.
        """
        pass

    @abstract_property
    def name(self):
        """Return method name used internally and in command line interface.

        The name should not include any spaces."""
        pass

    @abstractmethod
    def run(self, experiment: Experiment) -> MethodResult:
        pass
