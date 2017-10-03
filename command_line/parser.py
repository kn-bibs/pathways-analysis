import argparse
from collections import defaultdict
from copy import deepcopy


def group_arguments(args, group_names):
    """Group arguments into given groups + None group for all others"""
    groups = defaultdict(list)
    group = None
    for arg in args:
        if arg in group_names:
            group = arg
        else:
            groups[group].append(arg)
    return groups, groups[None]


class Argument:
    """Defines argument for `Parser`.

    In essence this is a wrapper for `argparse.add_argument()`,
    so most options (type, help) which work in standard python
    parser will work with Argument too. Additionally some nice
    feature, like automated naming are available.

    Worth to mention that when used with `MethodParser`,
    `type` and `help` will be automatically deduced.
    """

    def __init__(
            self, name=None, short=None, optional=True,
            as_many_as: 'Argument'=None, **kwargs
    ):
        """
        Args:
            name:
                overrides deduced argument name
            short:
                a single letter to be used as a short name
                (e.g. "c" will enable using "-c")
            optional:
                by default True, provide False
                to make the argument required
            as_many_as:
                if provided, will check if len() of the produced
                value is equal to len() of the provided argument
            **kwargs:
                other keyword arguments which are
                supported by `argparse.add_argument()`
        """
        self.name = name
        self.short_name = short
        self.optional = optional
        self.as_many_as = as_many_as
        self.kwargs = kwargs
        self.default = kwargs.get('default', None)

    @property
    def args(self):

        args = []
        if self.optional:
            if self.short_name:
                args.append(f'-{self.short_name}')
            args.append(f'--{self.name}')
        else:
            args.append(self.name)

        return args

    @staticmethod
    def as_numerous_as(myself, partner):
        # if we have callable, we can call it as many times as we want
        if callable(myself):
            return True
        if partner and myself:
            return len(partner) == len(myself)
        return True

    def validate(self, opts):
        myself = getattr(opts, self.name)

        if self.as_many_as:
            partner = getattr(opts, self.as_many_as.name)
            if not self.as_numerous_as(myself, partner):
                raise ValueError(
                    f'{self.name} for {len(myself)} {self.as_many_as.name} '
                    f'provided, expected for {len(partner)}'
                )


class Parser:
    """Parser is a wrapper around Python built-in `argparse.ArgumentParser`.

    Subclass the `Parser` to create your own parser.

    Use help, description and epilog properties to adjust the help screen.
    By default help and description will be auto-generated using docstring
    and defined arguments.

    Attach custom arguments and sub-parsers by defining class-variables
    with `Argument` and `Parser` instances.

    Example:

        class TheParser(Parser):
            help = 'This takes only one argument, but it is required'

            arg = Argument(optional=False, help='This is required')

        class MyParser(Parser):
            description = 'This should be a longer text'

            my_argument = Argument(type=int, help='some number')
            my_sub_parser = TheParser()

            epilog = 'You can create a footer with this'

        # To execute the parser use:

        parser = MyParser()

        # The commands will usually be `sys.argv[1:]`
        commands = '--my_argument 4 my_sub_parser value'.split()

        namespace = parser.parse(commands)

        # `namespace` is a normal `argparse.Namespace`
        assert namespace.my_argument == 4
        assert namespace.my_sub_parser.arg == 'value'

    Implementation details:

        To enable behaviour not possible with limited, plain `ArgumentParser`
        (e.g. to dynamically attach a sub-parser, or to chain two or more
        sub-parsers together) the stored actions and sub-parsers are:
            - not attached permanently to the parser,
            - attached in a tricky way to enable desired behaviour,
            - executed directly or in hierarchical order.

        Class-variables with parsers will be deep-copied on initialization,
        so you do not have to worry about re-use of parsers.
    """
    # sub-parsers will have dynamically populated name variable
    parser_name = None

    @property
    def help(self):
        """A short message, shown as summary on >parent< parser help screen.

        Help will be displayed for sub-parsers only.
        """
        return (
            ' Accepts: ' + ', '.join(self.arguments.keys())
        )

    @property
    def description(self):
        """Longer description of the parser.

        Description is shown when user narrows down the help
        to the parser with: `./run.py sub_parser_name -h`.
        """
        return (self.__doc__ or '').format(**vars(self))

    @property
    def epilog(self):
        """Use this to append text after the help message"""
        return None

    def __init__(self, parser_name=None, **kwargs):
        """Uses kwargs to pre-populate namespace of the `Parser`.

        Args:
            parser_name: a name used for identification of sub-parser
        """
        self.namespace = argparse.Namespace()
        self.parser_name = parser_name
        self.parser = argparse.ArgumentParser(description=self.description, epilog=self.epilog)

        self.arguments = {}
        # children parsers
        self.subparsers = {}

        # parses and arguments pulled up from children parsers
        self.lifted_parsers = {}
        self.lifted_args = {}

        attribute_handlers = {
            Argument: self.bind_argument,
            Parser: self.add_parser,
        }

        # register class attributes
        for name in dir(self):
            attribute = getattr(self, name)
            for attribute_type, handler in attribute_handlers.items():
                if isinstance(attribute, attribute_type):
                    handler(attribute, name)

        # initialize namespace
        for name, argument in self.all_arguments.items():
            setattr(self.namespace, name, argument.default)

        for name, value in kwargs.items():
            setattr(self.namespace, name, value)

        self.to_builtin_parser()
        self.kwargs = kwargs

    @property
    def all_subparsers(self):
        return {**self.subparsers, **self.lifted_parsers}

    @property
    def all_arguments(self):
        return {**self.arguments, **self.lifted_args}

    def to_builtin_parser(self):
        for argument in self.all_arguments.values():
            self.add_argument(argument)

    def add_argument(self, argument, parser=None):
        if not parser:
            parser = self.parser

        parser.add_argument(*argument.args, **argument.kwargs)

    def attach_subparsers(self):
        """Only in order to show a nice help, really

        There are some issues when using subparsers added with the built-in
        add_subparsers for parsing. Instead subparsers are handled in a
        custom implementation of parse_known_args (which really builds upon
        the built-in one, just tweaking some places).
        """
        # regenerate description and epilog: enables use of custom variables
        # (which may be not yet populated at init.) in descriptions epilogues
        self.parser.description = self.description
        self.parser.epilog = self.epilog

        native_sub_parser = self.parser.add_subparsers()

        for name, sub_parser in self.all_subparsers.items():

            if sub_parser.pull_to_namespace_above:
                continue

            parser = native_sub_parser.add_parser(
                help=sub_parser.help, name=name,
                description=sub_parser.description
            )

            for argument in sub_parser.arguments.values():
                self.add_argument(argument, parser)

    def add_parser(self, parser, name):
        # copy is needed as we do not want to share the class-method
        # stored data across all instances
        parser = deepcopy(parser)
        # and we want to use only this instance of the parser
        setattr(self, name, parser)
        parser.parser_name = name
        self.subparsers[name] = parser
        if parser.pull_to_namespace_above:
            self.lifted_args.update(parser.arguments)
            self.lifted_parsers.update(parser.subparsers)

    def __deepcopy__(self, memodict={}):
        return self.__class__(**self.kwargs)

    def bind_argument(self, argument, name=None):
        if not argument.name and name:
            argument.name = name
        self.arguments[name] = argument

    def parse_known_args(self, args):
        grouped_args, ungrouped_args = group_arguments(args, self.all_subparsers)

        for name, parser in self.subparsers.items():

            if parser.pull_to_namespace_above:

                namespace, not_parsed_args = parser.parse_known_args([
                    arg_str
                    for key in parser.subparsers
                    for arg_str in [key, *grouped_args[key]]
                    # only include the sub-parser if it was explicitly enlisted
                    if key in grouped_args
                ])

                for key, value in vars(namespace).items():
                    setattr(self.namespace, key, value)
            else:
                # only invoke sub-parser parsing if it was explicitly enlisted
                if name in grouped_args:
                    namespace, not_parsed_args = parser.parse_known_args(grouped_args[name])
                    setattr(self.namespace, name, namespace)
                else:
                    setattr(self.namespace, name, None)
                    not_parsed_args = None

            assert not not_parsed_args

        namespace, unknown_args = self.parser.parse_known_args(
            ungrouped_args,
            namespace=self.namespace
        )
        assert namespace is self.namespace

        # TODO: make validation (and production) errors raise in parser context
        self.validate(self.namespace)

        opts = self.produce(unknown_args)
        assert opts is self.namespace

        return self.namespace, unknown_args

    def validate(self, opts):
        if not opts:
            opts = self.namespace
        for argument in self.all_arguments.values():
            argument.validate(opts)

    @property
    def pull_to_namespace_above(self):
        return False

    def produce(self, unknown_args):
        """Post-process already parsed namespace.

        You can override this method to create a custom objects
        in the parsed namespace (e.g. if you cannot specify the
        target class with Argument(type=X), because X depends
        on two or more arguments).

        You can chery-pick the arguments which were not parsed
        by the current parser (e.g. when some step of parsing
        depends on provided arguments), but please remember
        to remove those from `unknown_args` list.

        Remember to operate on the provided list object (do not
        rebind the name with `unknown_args = []`, as doing so
        will have no effect: use `unknown_args.remove()` instead).
        """
        return self.namespace

    def parse(self, args):

        # Use the built-in help (just attach sub-parsers before).
        if '-h' in args or '--help' in args or not args:
            self.attach_subparsers()
            self.parser.parse_args(args)

        # Parse wisely, we need to support chaining sub-parsers,
        # validation and so on. Everything in parse_known_args.
        options, unknown_args = self.parse_known_args(args)
        assert not unknown_args

        return options

