import re
import subprocess as sp
import sys
import time

from .constants import MINIMUM_REQUIREMENTS

TECHNOLOGY_PARSER = re.compile(r'^(?P<name>\S+)')
VERSION_PARSERS = {
    requirement:
    re.compile(r'^{} ([0-9]+).([0-9]+).([0-9]+)'.format(requirement))
    for requirement in MINIMUM_REQUIREMENTS
}


class UnmetDependencyException(Exception):
    pass


def run_executable(
        command,
        stdin=None,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        wait=True,
        stream=False,
        quiet=False,
        returncode=0
):
    """Run a single shell command and wait for it to terminate.
    """
    command = [str(c) for c in command]
    if not quiet:
        print(' '.join(command), file=sys.stderr)
    p = sp.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        universal_newlines=wait
    )

    # Wait if desired.
    if wait:
        while p.poll() is None:
            if stream:
                for line in p.stdout:
                    print(line, end='')
                for line in p.stdout:
                    print(line, end='', file=sys.stderr)
            else:
                time.sleep(1)

        if not quiet and p.returncode != returncode:
            raise sp.CalledProcessError(p.returncode, ' '.join(command))

    return p


def run_chain(*commands, stdin=None, stdout=sp.PIPE, wait=True, stream=False):
    """Chain multiple commands by piping outputs to inputs.
    """
    assert len(commands) > 1
    processes = []
    for command in commands:
        _stdin = stdin
        _stdout = stdout
        _wait = wait
        _stream = stream
        if processes:
            _stdin = processes[-1].stdout
        if len(processes) != len(commands) - 1:
            _stdout = sp.PIPE
            _wait = False
            _stream = False
        p = run_executable(
            command, stdin=_stdin, stdout=_stdout, wait=_wait, stream=_stream
        )
        processes.append(p)

        if _stdin:
            _stdin.close()

    if wait:
        for p in processes:
            while p.poll() is None:
                time.sleep(1)
            if p.returncode != 0:
                raise sp.CalledProcessError(p.returncode, ' '.join(command))
    return processes


def get_version(requirement):
    p = run_executable([requirement], quiet=True, returncode=1)
    match = VERSION_PARSERS[requirement].match(p.stdout.read())
    return tuple(int(ver) for ver in match.groups()) if match else None


def check_dependencies():
    """Checks if executable dependencies have been met.
    """
    for requirement, minimum_version in MINIMUM_REQUIREMENTS.items():
        version = get_version(requirement)
        if version < minimum_version:
            raise UnmetDependencyException(
                '{} version {} is less than the minimum requirement {}'.
                format(requirement, version, minimum_version)
            )


def parse_technologies(lines):
    parsing = False
    technologies = set()
    for line in lines:
        if line.startswith('-'):
            parsing = True
            continue

        if parsing:
            if line.isspace():
                break
            match = TECHNOLOGY_PARSER.match(line)
            if match:
                technologies.add(match['name'])
    return technologies


def get_supported_technologies():
    """Runs 'kallisto bus --list' to fetch a list of supported technologies.
    """
    p = run_executable(['kallisto', 'bus', '--list'], quiet=True, returncode=1)
    return parse_technologies(p.stdout)


def download_whitelist(technology):
    """Downloads public whitelist barcodes for specified technology.
    """
    pass


def create_transcript_list(gtf_path, use_name=True, use_version=False):
    r = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if len(line) == 0 or line[0] == '#':
                continue
            l = line.strip().split('\t')  # noqa
            if l[2] == 'transcript':
                info = l[8]
                d = {}
                for x in info.split('; '):
                    x = x.strip()
                    p = x.find(' ')
                    if p == -1:
                        continue
                    k = x[:p]
                    p = x.find('"', p)
                    p2 = x.find('"', p + 1)
                    v = x[p + 1:p2]
                    d[k] = v

                if 'transcript_id' not in d or 'gene_id' not in d:
                    continue

                tid = d['transcript_id'].split(".")[0]
                gid = d['gene_id'].split(".")[0]
                if use_version:
                    if 'transcript_version' not in d or 'gene_version' not in d:
                        continue

                    tid += '.' + d['transcript_version']
                    gid += '.' + d['gene_version']
                gname = None
                if use_name:
                    if 'gene_name' not in d:
                        continue
                    gname = d['gene_name']

                if tid in r:
                    continue

                r[tid] = (gid, gname)
    return r
