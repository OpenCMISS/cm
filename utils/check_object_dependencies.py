#!/usr/bin/env python

"""
Compare object dependencies with the modules used

Must be run from the cm directory
"""

import os
import sys
import re


source_to_module = {}
module_to_source = {}
source_to_object = {}
modules_from_source = []

# Read source files and get module names and modules used
fortran_files = [f for f in os.listdir('src') if f.endswith('.f90')]
source_modules = {}
for fname in fortran_files:
    module_name_found = False
    used_modules = []
    for line in open('src' + os.sep + fname, 'r').readlines():
        if not module_name_found:
            match = re.match(r'MODULE\s+([a-zA-Z_]+)\s*', line, re.IGNORECASE)
            if match:
                module_name = match.group(1)
                module_name_found = True
        match = re.match(r'\s+USE\s+([a-zA-Z_]+)\s*', line, re.IGNORECASE)
        if match:
            used_modules.append(match.group(1).lower())
        if line.strip().upper() == 'CONTAINS':
            break
    if module_name_found:
        source_modules[fname] = used_modules
        source_to_module[fname] = module_name.lower()
        module_to_source[module_name.lower()] = fname
        modules_from_source.append(module_name.lower())
    else:
        sys.stderr.write("Couldn't find module name for %s\n" % fname)
        exit(1)

# Read object dependencies from Makefile
object_target_re = re.compile(r'\$\(OBJECT_DIR\)/([a-zA-Z_]+)\.o\s*:\s*'
                r'\$\(SOURCE_DIR\)/([a-zA-Z_]+\.[a-zA-Z0-9_]+).*')
makefile = open('Makefile', 'r').readlines()
objects = {}
object_makefile_contents = {}
state = 'searching'
for line in makefile:
    if state == 'searching':
        match = object_target_re.match(line)
        if match:
            object_name = match.group(1)
            if line.rstrip().endswith('\\'):
                state = 'in_object'
            else:
                object_makefile_contents[object_name] = []
            dependencies = []
            makefile_contents = []
            objects[object_name] = (match.group(2), dependencies)
            source_to_object[match.group(2)] = object_name + '.o'
    elif state == 'in_object':
        makefile_contents.append(line)
        match = re.match(r'\s*\$\(OBJECT_DIR\)/([a-zA-Z_]+)\.o\s*(\\)?', line)
        if match:
            dependencies.append(match.group(1))
        elif line.strip().startswith('$(FIELDML_OBJECT)'):
            dependencies.extend(['fieldml_util_routines',
                'fieldml_input_routines',
                'fieldml_output_routines',
                'fieldml_types'])
        elif line.strip().startswith('$(MACHINE_OBJECTS)'):
            dependencies.extend(['machine_constants_linux',
                    'machine_constants_win32',
                    'machine_constants_aix'])
        if not line.strip().endswith('\\'):
            state = 'searching'
            object_makefile_contents[object_name] = makefile_contents

# Compare the dependencies in the Makefile with
# the modules from the source file
missing_dependencies = {}
extra_dependencies = {}
for object in objects.keys():
    (source_file, dependencies) = objects[object]
    print('* %s' % source_file)
    if not source_file.endswith('.f90'):
        print("  Skipping non Fortran file.")
        continue
    try:
        # Skip external modules (not in modules_from_source)
        modules_used = set([
            m.lower() for m in source_modules[source_file]
            if m.lower() in modules_from_source])
    except KeyError:
        print("  Didn't find source file.")
        continue
    dependency_sources = []
    for d in dependencies:
        try:
            dependency_sources.append(objects[d][0])
        except KeyError:
            # Dependency not listed as a target itself
            if os.path.isfile('src/' + d + '.f90'):
                dependency_sources.append(d + '.f90')
            else:
                print("  Skipping dependency %s" % d)
    # convert the list of source file dependencies to module names,
    # skipping non-fortran files
    dependency_modules = set([
        source_to_module[d].lower() for d in dependency_sources
        if d.endswith('.f90')])
    for d in dependency_sources:
        if not d.endswith('.f90'):
            print("  Can't check dependency %s" % d)
    not_in_source = dependency_modules - modules_used
    not_in_makefile = modules_used - dependency_modules
    extra_dependencies[object] = not_in_source
    missing_dependencies[object] = not_in_makefile
    if not_in_source:
        print("  Makefile dependency not used:")
        print("    - " + "\n    - ".join(not_in_source))
    if not_in_makefile:
        print("  Modules used but not listed as a dependency:")
        print("    + " + "\n    + ".join(not_in_makefile))


def update_contents(obj, contents, missing, extra):
    new_content = []
    for content in contents:
        content = content.strip()
        is_extra = False
        for extra_module in extra:
            extra_obj = source_to_object[module_to_source[extra_module]]
            if content.find('/' + extra_obj + '.') > -1:
                is_extra = True
        if not is_extra:
            new_content.append(content)
    for missing_module in missing:
        missing_obj = source_to_object[module_to_source[missing_module]]
        new_content.append("$(OBJECT_DIR)/%s" % missing_obj)
    new_content = ['\t' + l.rstrip('\\').strip() + ' \\' for l in new_content]
    new_content.sort()
    try:
        new_content[-1] = '\t' + new_content[-1].rstrip('\\').strip()
    except IndexError:
        pass
    return new_content


if "-fix" in sys.argv:
    # Make an attempt to fix problems, often doesn't work
    for object in object_makefile_contents.keys():
        try:
            object_makefile_contents[object] = update_contents(
                    object,
                    object_makefile_contents[object],
                    missing_dependencies[object],
                    extra_dependencies[object])
        except KeyError:
            object_makefile_contents[object] = [
                    l.rstrip() for l in object_makefile_contents[object]]
            # Non Fortran file or other skipped object
            continue
    orig_makefile = makefile
    # Read object dependencies from Makefile
    makefile = open('Makefile', 'w')
    state = 'searching'
    for line in orig_makefile:
        if state == 'searching':
            makefile.write(line)
            match = object_target_re.match(line)
            if match:
                state = 'in_object'
                object_name = match.group(1)
                makefile.write('\n'.join(
                    object_makefile_contents[object_name]))
                makefile.write('\n')
        elif state == 'in_object':
            if not line.strip().endswith('\\'):
                state = 'searching'
    makefile.close()
