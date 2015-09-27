#!/usr/bin/env python
"""
Update OpenCMISS Fortran code to use the new cmfe_ prefix.
Also use an array to set the equation set and problem specifications.
"""

import os
import re
import sys


# Pairs of function and regex to search for
converters = []
readers = []


def converter(regex):
    """Register function that converts a regex match
    """
    def decorator(fn):
        converters.append((fn, re.compile(regex, re.IGNORECASE)))
        return fn
    return decorator


def reader(regex):
    """Register function that reads information from a regex match
    """
    def decorator(fn):
        readers.append((fn, re.compile(regex, re.IGNORECASE)))
        return fn
    return decorator


@reader(r'CALL\s+CMISSProblem_SpecificationSet\((.*)\)')
def read_problem_specification(match, info):
    if 'problem_specification' in info:
        raise ValueError("Multiple problem specifications found")
    args = [a.strip() for a in match.group(1).split(',')]
    info['problem_specification'] = args[1:4]


@converter(r'CALL\s+CMISSProblem_CreateStart\((.*)\)')
def convert_problem_start(match, info):
    args = [a.strip() for a in match.group(1).split(',')]
    problem_spec = '[%s]' % ','.join(info['problem_specification'])
    new_args = args[:1] + [problem_spec] + args[1:]
    return "CALL CMISSProblem_CreateStart(%s)" % ','.join(new_args)


@converter(r'CALL\s+CMISSEquationsSet_CreateStart\((.*)\)')
def convert_equations_set_start(match, info):
    args = [a.strip() for a in match.group(1).split(',')]
    spec = '[%s]' % ','.join(args[3:6])
    new_args = args[:3] + [spec] + args[6:]
    return "CALL CMISSEquationsSet_CreateStart(%s)" % ','.join(new_args)


@converter(r'CALL\s+CMISSProblem_SpecificationSet\(.*\)')
def remove_problem_specification(match, info):
    return ""


@converter(r'!Set the problem to be.*$')
def remove_problem_spec_comment(match, info):
    return ""


def convert(source):
    """Convert source file, generating new lines
    """
    # Read all of input as we need to go over it twice
    lines = source.readlines()

    # Find information
    info = {}
    for line, _ in full_lines(lines):
        for read, regex in readers:
            match = regex.search(line)
            if match:
                read(match, info)

    # Now convert file
    for line, orig_lines in full_lines(lines):
        for convert, regex in converters:
            match = regex.search(line)
            if match:
                converted = convert(match, info)
                yield fix_lines(
                        line[:match.start()] + converted + line[match.end():])
                break
        else:
            yield orig_lines


def full_lines(source):
    """ Read full lines over line continuations,
        yielding the parsed line and the original content
    """
    full_line = ""
    original_lines = ""
    for line in source:
        original_lines += line
        full_line += line.replace('&', '').replace('\n', '')
        if not line.rstrip().endswith('&'):
            yield (full_line, original_lines)
            full_line = ""
            original_lines = ""


def get_indent(line):
    indent = line.replace(line.lstrip(), '')
    return indent


def fix_lines(lines):
    if lines.strip() == "":
        return ""
    return '\n'.join(fix_line(line) for line in lines.splitlines()) + '\n'


def fix_line(line):
    """ Add Fortran line continuation after a comma if line is too long
    """
    MAX_LEN = 132
    if len(line) <= MAX_LEN:
        return line

    indent = get_indent(line)
    remaining = line.rstrip()
    content = ""
    while len(remaining) > MAX_LEN:
        break_pos = remaining.rfind(',', 0, 130) + 1
        if break_pos < 0:
            raise ValueError("Couldn't split line: %s\n" % line)
        content += remaining[0:break_pos] + ' &\n'
        remaining = indent + '  & ' + remaining[break_pos:]
    content += remaining

    return content


constant_prefix_re = re.compile(r'\bCMISS_')
prefix_re = re.compile(r'\bCMISS([A-Za-z])')
type_suffix_re = re.compile(r'([0-9])_CMISS')


def convert_prefix(line):
    # Convert constants, subroutine calls and type suffixes
    line = constant_prefix_re.sub('CMFE_', line)
    line = prefix_re.sub(r'cmfe_\1', line)
    line = type_suffix_re.sub(r'\1_CMFE', line)
    # Fix up stuff like CMISSDP etc that shouldn't be changed
    line = line.replace('cmfe_DP', 'CMFEDP')
    line = line.replace('cmfe_SP', 'CMFESP')
    line = line.replace('cmfe_QP', 'CMFEQP')
    line = line.replace('cmfe_DPC', 'CMFEDPC')
    line = line.replace('cmfe_SPC', 'CMFESPC')
    line = line.replace('cmfe_Intg', 'CMFEIntg')
    line = line.replace('cmfe_SIntg', 'CMFESIntg')
    line = line.replace('cmfe_LIntg', 'CMFELIntg')
    return line


if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[1] + '.fixed'
    with open(infile, 'r') as input:
        with open(outfile, 'w') as output:
            for newline in convert(input):
                output.write(convert_prefix(newline))
    os.rename(outfile, infile)
