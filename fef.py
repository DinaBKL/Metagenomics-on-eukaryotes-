#!/usr/bin/env python3

import argparse
import glob
import subprocess
import sys
from datetime import datetime

USAGE = """
  * fef (Python version) *

  Construct a script using expressions.
  By default it just prints the script on the screen.

  OPTIONS:
    -c, --command     The command line template (supports @ for files, # for numbers)
    -e, --expression  Expression to expand (e.g., "*.fna")
    -o, --outfile     Save the result to a file
    -n, --number      Number of commands to generate if no expression is given
    -r, --run         Run the commands one by one (asks confirmation)
    --repeat          Repeat the generated commands N times

  Example:
    fef -c "formatdb -i @ -p F" -e "*.fna"
    fef -c "mkdir temp#" -n 5
"""


def expand_commands(command, expression=None, number=0):
    commands = []
    if expression:
        files = sorted(glob.glob(expression))
        for idx, f in enumerate(files, start=1):
            new_cmd = command.replace("@", f)
            new_cmd = new_cmd.replace("#", str(idx))
            commands.append(new_cmd)
    elif number > 0:
        for idx in range(1, number + 1):
            new_cmd = command.replace("#", str(idx))
            commands.append(new_cmd)
    return commands


def write_to_file(outfile, command, expression, repeat, commands):
    with open(outfile, "w") as out:
        date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        out.write(f"# SCRIPT GENERATED USING fef on {date}\n")
        if repeat:
            out.write(f"# command = {command} , expression = {expression} , repeat = {repeat}\n")
        else:
            out.write(f"# command = {command} , expression = {expression}\n")
        out.write(f"# A total of {len(commands)} commands are present in this file\n")
        for cmd in commands:
            out.write(cmd + "\n")


def run_commands(commands):
    input("Going to run the commands now, press ENTER to start or Control-C to quit")
    total = len(commands)
    for i, cmd in enumerate(commands, start=1):
        print(f"Running [{i} of {total}]: {cmd}")
        subprocess.run(cmd, shell=True)


def main():
    parser = argparse.ArgumentParser(usage=USAGE, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--command", required=True, help="Command template")
    parser.add_argument("-e", "--expression", help="File expression to expand")
    parser.add_argument("-n", "--number", type=int, default=0, help="Number of commands")
    parser.add_argument("--repeat", type=int, default=0, help="Repeat commands N times")
    parser.add_argument("-o", "--outfile", help="Output file to save commands")
    parser.add_argument("-r", "--run", action="store_true", help="Run commands instead of just printing")
    args = parser.parse_args()

    if not args.expression and args.number == 0:
        print(USAGE)
        sys.exit("Error: Must specify either -e (expression) or -n (number)")

    commands = expand_commands(args.command, args.expression, args.number)

    # Handle repetition
    if args.repeat > 0:
        commands = commands * args.repeat

    # Print c
