import argparse
from collections import defaultdict


def group_arguments(args, group_names):
    """Group arguments into given groups + None group for all others

    Args:
        args:
        group_names:

    Returns:

    """
    groups = defaultdict(list)
    group = None
    for arg in args:
        if arg in group_names:
            group = arg
        else:
            groups[group].append(arg)
    return groups, groups[None]


class Argument:

    def __init__(self, **kwargs):
        # store the name if provided
        self.name = kwargs.pop('name', None)
        self.short_name = kwargs.pop('short', None)
        # parse some custom kwargs and froze arguments. In signature? In partial?
        self.kwargs = kwargs

    @property
    def args(self):
        if self.short_name:
            return [f'-{self.short_name}', f'--{self.name}']
        return [f'--{self.name}']


class Parser:
    """
    Uses arguments and subparser defined as class properties to initialize an ArgumentParser.
    """
    # sub-parsers will have dynamically populated name variable
    name = None

    @property
    def help(self):
        return (
            self.__doc__.format(name=self.name) +
            ' Accepts: ' + ', '.join(self.arguments.keys())
        )

    @property
    def description(self):
        return self.__doc__.format(name=self.name)

    def __init__(self, **kwargs):
        """Uses kwargs to populate namespace of the Parser."""
        self.namespace = argparse.Namespace()
        self.parser = argparse.ArgumentParser()

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
        for name, attribute in vars(self.__class__).items():
            for attribute_type, handler in attribute_handlers.items():
                if isinstance(attribute, attribute_type):
                    handler(attribute, name)

        # initialize namespace
        for name in self.all_arguments:
            setattr(self.namespace, name, None)

        for name, value in kwargs.items():
            setattr(self.namespace, name, value)

        self.to_builtin_parser()

    @property
    def all_subparsers(self):
        return {**self.subparsers, **self.lifted_parsers}

    @property
    def all_arguments(self):
        return {**self.arguments, **self.lifted_args}

    def to_builtin_parser(self):
        native_sub_parser = self.parser.add_subparsers()

        for argument in self.all_arguments.values():
            self.parser.add_argument(*argument.args, **argument.kwargs)

        for name, sub_parser in self.all_subparsers.items():

            if sub_parser.pull_to_namespace_above:
                continue

            parser = native_sub_parser.add_parser(
                help=sub_parser.help, name=name,
                description=sub_parser.description
            )

            for argument in sub_parser.arguments.values():
                parser.add_argument(*argument.args, **argument.kwargs)

    def add_parser(self, parser, name):
        parser.name = name
        self.subparsers[name] = parser
        if parser.pull_to_namespace_above:
            self.lifted_args.update(parser.arguments)
            self.lifted_parsers.update(parser.subparsers)

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
                ])

                for key, value in vars(namespace).items():
                    setattr(self.namespace, key, value)
            else:
                namespace, not_parsed_args = parser.parse_known_args(grouped_args[name])
                setattr(self.namespace, name, namespace)

            assert not not_parsed_args

        namespace, unknown_args = self.parser.parse_known_args(
            ungrouped_args,
            namespace=self.namespace
        )
        assert namespace is self.namespace

        opts = self.produce(unknown_args)

        return self.namespace, unknown_args

    def parse(self, args):
        return self.parser.parse_args(args, namespace=self.namespace)


    @property
    def pull_to_namespace_above(self):
        """If you do not want a parser to require --name, set this to True"""
        return False

    def produce(self, unknown_args):
        """Hook to modify arguments after they were supplied"""
        return unknown_args


