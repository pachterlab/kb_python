import subprocess as sp
from unittest import mock, TestCase
from unittest.mock import call

import kbcount.utils as utils


class TestUtils(TestCase):
    def test_run_executable(self):
        p = utils.run_executable(['echo', 'TEST'])
        self.assertEqual(p.stdout.read(), 'TEST\n')

    def test_run_exectuable_raises_exception(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_executable(['bash', 'exit', '1'])

    def test_run_exectuable_with_returncode(self):
        utils.run_executable(['bash', 'exit', '1'], returncode=1)

    def test_run_executable_no_wait(self):
        with mock.patch('kbcount.utils.sp') as sp_mock:
            sp_mock.Popen().returncode = 0
            utils.run_executable(['echo', 'TEST'], wait=False)
            sp_mock.Popen().poll.assert_not_called()

    def test_run_executable_with_stream(self):
        with mock.patch('kbcount.utils.print') as print_mock:
            utils.run_executable(['echo', 'TEST'], stream=True)
            print_mock.assert_has_calls([
                call('echo TEST'),
                call('TEST\n', end=''),
            ])

    def test_run_chain(self):
        ps = utils.run_chain(['echo', 'TEST'], ['grep', 'T'])
        self.assertEqual(ps[1].stdout.read(), 'TEST\n')

    def test_run_chain_fails_single_command(self):
        with self.assertRaises(AssertionError):
            utils.run_chain(['echo', 'TEST'])

    def test_run_chain_raises_exception_when_dead(self):
        with self.assertRaises(sp.SubprocessError):
            utils.run_chain(['sleep', '5'], ['grep', 'TEST'], ['ls'])
