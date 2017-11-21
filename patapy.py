#!/usr/bin/python3.6
import sys
from command_line import CLI
from methods import Method, MethodResult


# TODO: HTML version with Bootstrap and bootstrap-table;
# f-strings should work as simple-templates; option "show in browser"?
# or we could use jupiter notebooks?
def render_text_table(method: Method, results: MethodResult):
    print(f'Results of {method.name} run')

    print(results.description)

    print('\t'.join(results.columns))

    for result in results.scored_list:
        for column in results.columns:
            v = getattr(result, column)
            print(v, end='\t')
        print()

    if results.files:
        print('There are additional output files in following locations:')
        print(results.files)

    # TODO: method.parameters?


def run(argv):
    args = CLI().parse_args(argv[1:])
    results = args.method.run(args.experiment)
    render_text_table(args.method, results)
    return results


if __name__ == '__main__':
    run(sys.argv)
