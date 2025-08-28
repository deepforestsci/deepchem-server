import os
import unittest

import pytest

from deepchem_server.core.address import DeepchemAddress


class TestDeepchemAddress(unittest.TestCase):
    """
    Tests for DeepchemAddress

    """

    def setUp(self):
        raw_address1 = "namespace/username/filename"
        raw_address2 = "deepchem://namespace/username/filename"
        raw_address3 = "deepchem://namespace/username/nested/filename"
        raw_address4 = "deepchem://namespace/filename"
        self.address1 = DeepchemAddress(raw_address1)
        self.address2 = DeepchemAddress(raw_address2, 'data')
        self.address3 = DeepchemAddress(raw_address3, 'data')
        with pytest.raises(ValueError) as excinfo:
            DeepchemAddress(raw_address4, 'data')
        self.assertEqual(str(excinfo.value), 'Invalid deepchem address')

    def test_get_key(self):
        self.assertEqual(DeepchemAddress.get_key(self.address1.address),
                         'filename')
        self.assertEqual(DeepchemAddress.get_key(self.address2.address),
                         'filename')
        self.assertEqual(DeepchemAddress.get_key(self.address3.address),
                         'nested/filename')

    def test_parse(self):
        self.assertEqual(DeepchemAddress.parse_address(self.address1.address), {
            'profile': 'namespace',
            'project': 'username',
            'key': 'filename'
        })
        self.assertEqual(DeepchemAddress.parse_address(self.address2.address), {
            'profile': 'namespace',
            'project': 'username',
            'key': 'filename'
        })
        self.assertEqual(DeepchemAddress.parse_address(self.address3.address), {
            'profile': 'namespace',
            'project': 'username',
            'key': 'nested/filename'
        })

    def test_get_path(self):
        parsed_address1 = DeepchemAddress.parse_address(self.address1.address)
        parsed_address2 = DeepchemAddress.parse_address(self.address2.address)
        parsed_address3 = DeepchemAddress.parse_address(self.address3.address)

        local_format1 = os.path.join(parsed_address1['profile'],
                                     parsed_address1['project'],
                                     parsed_address1['key'])
        local_format2 = os.path.join(parsed_address2['profile'],
                                     parsed_address2['project'],
                                     parsed_address2['key'])
        local_format3 = os.path.join(parsed_address3['profile'],
                                     parsed_address3['project'],
                                     parsed_address3['key'])

        local_storage_loc1 = os.path.join(parsed_address1['profile'],
                                          parsed_address1['project'])
        local_storage_loc2 = os.path.join(parsed_address2['profile'],
                                          parsed_address2['project'])
        local_storage_loc3 = os.path.join(parsed_address3['profile'],
                                          parsed_address3['project'])

        # local format
        self.assertEqual(
            DeepchemAddress.get_path(local_storage_loc1,
                                     self.address1.address,
                                     format='local'), local_format1)
        self.assertEqual(
            DeepchemAddress.get_path(local_storage_loc2,
                                     self.address2.address,
                                     format='local'), local_format2)
        self.assertEqual(
            DeepchemAddress.get_path(local_storage_loc3,
                                     self.address3.address,
                                     format='local'), local_format3)

        # local format with only key
        self.assertEqual(
            DeepchemAddress.get_path(local_storage_loc1,
                                     DeepchemAddress.get_key(
                                         self.address1.address),
                                     format='local'), local_format1)
        self.assertEqual(
            DeepchemAddress.get_path(local_storage_loc2,
                                     DeepchemAddress.get_key(
                                         self.address2.address),
                                     format='local'), local_format2)
        self.assertEqual(
            DeepchemAddress.get_path(local_storage_loc3,
                                     DeepchemAddress.get_key(
                                         self.address3.address),
                                     format='local'), local_format3)

    def test_get_parent_key(self):
        self.assertEqual(DeepchemAddress.get_parent_key(self.address1.address),
                         '')
        self.assertEqual(DeepchemAddress.get_parent_key(self.address2.address),
                         '')
        self.assertEqual(DeepchemAddress.get_parent_key(self.address3.address),
                         'nested/')

    def test_get_object_name(self):
        self.assertEqual(DeepchemAddress.get_object_name(self.address1.address),
                         'filename')
        self.assertEqual(DeepchemAddress.get_object_name(self.address2.address),
                         'filename')
        self.assertEqual(DeepchemAddress.get_object_name(self.address3.address),
                         'filename')

    def test_repr(self):
        self.assertEqual(repr(self.address1),
                         'deepchem://namespace/username/filename')
        self.assertEqual(repr(self.address2),
                         'deepchem://namespace/username/filename')
