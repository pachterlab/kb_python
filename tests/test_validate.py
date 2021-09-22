import os
from unittest import mock, TestCase
from unittest.mock import call

import kb_python.validate as validate
from tests.mixins import TestMixin


@validate.validate_files()
def dummy(*args, **kwargs):
    return


@validate.validate_files()
def dummy_str(*args, **kwargs):
    return 'test'


@validate.validate_files()
def dummy_dict(*args, **kwargs):
    return {'test': 'testfile'}


@validate.validate_files()
def dummy_tuple(*args, **kwargs):
    return 'test1', 'test2'


class TestValidate(TestMixin, TestCase):

    def test_validate_bus(self):
        validate.validate_bus(self.bus_path)

    def test_validate_bus_failed_parse(self):
        with mock.patch('kb_python.validate.run_executable') as run_executable:
            run_executable.return_value = None, '', None
            with self.assertRaises(validate.ValidateError):
                validate.validate_bus('path')

    def test_validate_bus_no_records(self):
        with mock.patch('kb_python.validate.run_executable') as run_executable:
            run_executable.return_value = None, 'Read in 0 BUS records', None
            with self.assertRaises(validate.ValidateError):
                validate.validate_bus('path')

    def test_validate_mtx(self):
        validate.validate_mtx(self.matrix_path)

    def test_validate_mtx_raises_on_error(self):
        with mock.patch('kb_python.validate.scipy.io.mmread') as mmread:
            mmread.side_effect = ValueError('test')
            with self.assertRaises(validate.ValidateError):
                validate.validate_mtx('path')

    def test_validate(self):
        mock_validators = {
            '.bus': mock.MagicMock(),
            '.mtx': mock.MagicMock(),
        }
        bus_path = os.path.join(self.temp_dir, 'bus.bus')
        open(bus_path, 'w')
        mtx_path = os.path.join(self.temp_dir, 'mtx.mtx')
        open(mtx_path, 'w')

        with mock.patch('kb_python.validate.VALIDATORS', mock_validators):
            validate.validate(bus_path)
            mock_validators['.bus'].assert_called_once_with(bus_path)
            validate.validate(mtx_path)
            mock_validators['.mtx'].assert_called_once_with(mtx_path)

    def test_validate_doesnt_exist(self):
        with self.assertRaises(validate.ValidateError):
            validate.validate('nonexistent/path')

    def validate_files(self):
        with mock.patch('kb_python.validate.validate') as v,\
            mock.patch('kb_python.validate.os.path.exists') as exists:
            exists.return_value = True
            self.assertIsNone(dummy('f1', 1, kwarg1='f2', kwarg2=2))
            self.assertEqual(2, v.call_count)
            v.assert_has_calls([call('f1'), call('f2')])

    def validate_files_str(self):
        with mock.patch('kb_python.validate.validate') as v,\
            mock.patch('kb_python.validate.os.path.exists') as exists:
            exists.return_value = True
            self.assertEqual('test', dummy_str('f1', 1, kwarg1='f2', kwarg2=2))
            self.assertEqual(3, v.call_count)
            v.assert_has_calls([call('f1'), call('f2'), call('test')])

    def validate_files_dict(self):
        with mock.patch('kb_python.validate.validate') as v,\
            mock.patch('kb_python.validate.os.path.exists') as exists:
            exists.return_value = True
            self.assertEqual({'test': 'testfile'},
                             dummy_str('f1', 1, kwarg1='f2', kwarg2=2))
            self.assertEqual(3, v.call_count)
            v.assert_has_calls([call('f1'), call('f2'), call('testfile')])

    def validate_files_tuple(self):
        with mock.patch('kb_python.validate.validate') as v,\
            mock.patch('kb_python.validate.os.path.exists') as exists:
            exists.return_value = True
            self.assertEqual(('test1', 'test2'),
                             dummy_str('f1', 1, kwarg1='f2', kwarg2=2))
            self.assertEqual(4, v.call_count)
            v.assert_has_calls([
                call('f1'),
                call('f2'),
                call('test1'),
                call('test2')
            ])

    def test_validate_off(self):
        with mock.patch('kb_python.validate.is_validate') as is_validate,\
            mock.patch('kb_python.validate.validate_bus') as v:
            is_validate.return_value = False
            validate.validate('path/to/bus.bus')
            v.assert_not_called()
