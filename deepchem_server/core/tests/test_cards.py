import unittest
import tempfile
import os
from deepchem_server.core.cards import DataCard, ModelCard


class TestCards(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.tempfile_path = os.path.join(self.tempdir.name, 'temp.txt')
        self.address = "deepchem://namespace/username/working_dir/filename.extension"
        self.description = "this is a random description"
        self.kwargs = {'one': 'one', 'two': 'two'}

    def test_data_card(self):
        card = DataCard(address=self.address,
                        file_type="csv",
                        data_type="pandas.DataFrame",
                        feat_kwargs=self.kwargs,
                        other_meta_data="metadata",
                        shape=[10, 5])

        assert isinstance(card.shape, list)
        card_as_json = card.to_json()
        with open(self.tempfile_path, 'w') as fp:
            fp.write(card_as_json)

        with open(self.tempfile_path, 'r') as fp:
            card_as_json_from_file = fp.readlines()

        read_card = DataCard.from_json(card_as_json_from_file[0])
        for key, value in card.__dict__.items():
            assert value == getattr(read_card, key)

    def test_model_card(self):
        model_type = 'random_forest_classifier'
        train_dataset_address = "deepchem://namespace/username/working_dir/train_dataset.address"
        card = ModelCard(address=self.address,
                         model_type=model_type,
                         init_kwargs=self.kwargs,
                         train_dataset_address=train_dataset_address)
        with open(self.tempfile_path, 'w') as fp:
            fp.write(card.to_json())

        with open(self.tempfile_path, 'r') as fp:
            card_as_json_from_file = fp.readlines()

        read_card = ModelCard.from_json(card_as_json_from_file[0])
        for key, value in card.__dict__.items():
            assert value == getattr(read_card, key)
