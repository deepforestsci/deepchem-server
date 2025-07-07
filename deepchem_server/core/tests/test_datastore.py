import pytest
import os
import tempfile
import pandas as pd
import numpy as np
import deepchem as dc
from PIL import Image
from PIL.PngImagePlugin import PngImageFile
from deepchem_server.core.datastore import DiskDataStore
from deepchem_server.core import cards
from deepchem_server.core.cards import DataCard, ModelCard


def test_disk_datastore_in_memory_upload(disk_datastore):
    """Test upload in memory."""
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')
    df = pd.DataFrame({"foo": [1, 2], "bar": [3, 4]})
    data_address = disk_datastore.upload_data_from_memory(
        df, "test.csv", data_card)
    df2 = disk_datastore.get(data_address)
    assert len(df2) == 2
    assert data_address == "deepchem://test/user/test.csv"


def test_disk_datastore_in_memory_nested_upload(disk_datastore):
    """Test uploading to a nested folder in memory"""
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')
    df = pd.DataFrame({"foo": [1, 2], "bar": [3, 4]})
    data_address = disk_datastore.upload_data_from_memory(
        df, "test upload from memory/test.csv", data_card)
    df2 = disk_datastore.get(data_address)
    assert len(df2) == 2
    assert data_address == "deepchem://test/user/test upload from memory/test.csv"


def test_diskdatastore_in_memory_disk_dataset_upload(disk_datastore):
    """Test uploading and downloading dataset."""
    X = np.random.rand(10, 10)
    y = np.random.rand(10)
    data_card = cards.DataCard(address='',
                               file_type='dir',
                               data_type='dc.data.DiskDataset',
                               description='this is a disk dataset')
    data = dc.data.NumpyDataset(X, y)
    data_address = disk_datastore.upload_data_from_memory(
        data, 'dataset', data_card)
    data2 = disk_datastore.get(data_address)
    assert isinstance(data2, dc.data.DiskDataset)
    assert np.isclose(X, data2.X).all()


def test_datastore_prepopulated(tmpdir):
    """Test ability to connect datastore to existing directory."""
    base = str(tmpdir)
    session_dir = os.path.join(base, 'profile', 'project')
    os.makedirs(session_dir)
    df = pd.DataFrame({'foo': [1, 2], 'bar': [3, 4]})
    datapath = os.path.join(session_dir, 'test.csv')
    df.to_csv(datapath)
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a csv file')
    data_card_path = os.path.join(session_dir, 'test.csv.cdc')
    with open(data_card_path, 'wb') as f:
        f.write(bytes(data_card))

    dfs = DiskDataStore(profile_name='profile',
                        project_name='project',
                        basedir=base)
    test_address = "deepchem://profile/project/test.csv"
    df_get = dfs.get(test_address)  # noqa


def test_disk_datastore_download_object(disk_datastore, tmp_csv, tmpdir):
    # test file download
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')

    address = disk_datastore.upload_data('tmpcsv.csv', tmp_csv, data_card)

    tmpfilename = os.path.join(tmpdir, 'temp.csv')
    disk_datastore.download_object(address, tmpfilename)
    assert os.path.isfile(tmpfilename)

    # test dir download
    disk_datastore.add_dir('tmp_folder')
    disk_datastore.upload_data('tmp_folder/tmpcsv.csv', tmp_csv, data_card)
    folder_path = os.path.join(tmpdir, 'tmp_folder1')
    os.mkdir(folder_path)
    address = disk_datastore.storage_loc + '/' + 'tmp_folder'
    disk_datastore.download_object(address, folder_path)
    assert os.path.isdir(folder_path)
    assert os.path.isfile(os.path.join(folder_path, 'tmpcsv.csv'))


def test_disk_datastore_upload_data_from_memory(disk_datastore):
    """Test uploading data."""
    X = np.random.rand(10, 10)
    y = np.random.rand(10)
    data_card = cards.DataCard(address='',
                               file_type='dir',
                               data_type='dc.data.DiskDataset',
                               description='this is a disk dataset')
    data = dc.data.NumpyDataset(X, y)
    data_address = disk_datastore.upload_data_from_memory(
        data, 'dataset', data_card)
    data2 = disk_datastore.get_data(data_address)
    assert isinstance(data2, dc.data.DiskDataset)
    np.testing.assert_allclose(X, data2.X)


def test_disk_datastore_get_csv_data(disk_datastore, tmp_csv):
    """Test getting csv file_type data."""
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')

    address = disk_datastore.upload_data('tmpcsv.csv', tmp_csv, data_card)
    assert address == 'deepchem://test/user/tmpcsv.csv'

    df_retrieved = disk_datastore.get_data(address)
    assert df_retrieved.equals(pd.read_csv(tmp_csv))

    card_retrieved = disk_datastore.get(address + '.cdc')
    assert card_retrieved.address == data_card.address
    assert card_retrieved.file_type == 'csv'
    assert card_retrieved.data_type == 'pandas.DataFrame'


def test_disk_datastore_get_sample_csv_data(disk_datastore, tmp_csv):
    """Test getting csv file_type data sample."""
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')

    address = disk_datastore.upload_data('tmpcsv.csv', tmp_csv, data_card)
    assert address == 'deepchem://test/user/tmpcsv.csv'

    assert hasattr(disk_datastore, 'sample_rows')
    assert disk_datastore.sample_rows == 100

    # for testing
    disk_datastore.sample_rows = 2

    df_retrieved = disk_datastore.get(address, fetch_sample=True)
    assert df_retrieved.equals(pd.read_csv(tmp_csv, nrows=2))

    card_retrieved = disk_datastore.get(address + '.cdc', fetch_sample=True)
    assert card_retrieved.address == data_card.address
    assert card_retrieved.file_type == 'csv'
    assert card_retrieved.data_type == 'pandas.DataFrame'


def test_disk_datastore_get(disk_datastore):
    """Test getting card."""
    data_card = cards.DataCard(address='', file_type='json', data_type='json')
    address = disk_datastore.upload_data_from_memory(data_card,
                                                     'tempcard.json.cdc',
                                                     data_card).strip("'")
    assert address == 'deepchem://test/user/tempcard.json.cdc'

    card_retrieved = disk_datastore.get(address)
    assert card_retrieved.address == data_card.address


def test_png_upload_get(disk_datastore):
    """
    Test upload and get of png file in disk datastore
    """
    tempdir = tempfile.TemporaryDirectory()
    base_path = os.path.dirname(os.path.abspath(__file__))
    im = Image.open(os.path.join(base_path, 'assets/test-image.png'))
    path = os.path.join(tempdir.name, 'file.png')
    im.save(path)
    card = DataCard(address='', file_type='png', data_type='png')
    image_address = disk_datastore.upload_data_from_memory(
        im, 'temp_image.png', card)

    image_returned = disk_datastore.get(image_address)
    assert isinstance(image_returned, PngImageFile)
    assert (np.array(image_returned) == np.array(im)).all()


def test_txt_upload_get_delete(disk_datastore):
    """
    Test upload and get of txt file in disk datastore
    """
    # The file is present in the assets folder
    base_path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_path, 'assets/log.txt')
    card = DataCard(address='', file_type='txt', data_type='text/plain')
    file_address = disk_datastore.upload_data('log.txt', path, card)

    file_address = disk_datastore.upload_data('log.txt', path, card)

    file_returned = disk_datastore.get_data(file_address)
    assert (file_returned is not None)
    assert (isinstance(file_returned, list))
    with open(path, 'r') as f:
        assert (f.readlines() == file_returned)

    disk_datastore.delete_object(file_address, kind='data')

    with pytest.raises(FileNotFoundError):
        file_returned = disk_datastore.get_data(file_address)


def test_xml_upload_get_delete(disk_datastore):
    """
    Test upload and get of XML file in disk datastore
    """
    # The file is present in the assets folder
    base_path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_path, 'assets/sample.xml')
    card = DataCard(address='', file_type='xml', data_type='text/plain')
    file_address = disk_datastore.upload_data('sample.xml', path, card)

    file_returned = disk_datastore.get_data(file_address)
    assert (file_returned is not None)
    assert (isinstance(file_returned, list))
    with open(path, 'r') as f:
        assert (f.readlines() == file_returned)

    disk_datastore.delete_object(file_address, kind='data')

    with pytest.raises(FileNotFoundError):
        file_returned = disk_datastore.get_data(file_address)


def test_get_file_size_disk_datastore(disk_datastore):
    """
    Test get_file_size_in_bytes method
    """
    tempdir = tempfile.TemporaryDirectory()
    df = pd.DataFrame({'test': [0, 1, 2, 3, 4]})
    path = os.path.join(tempdir.name, 'file.csv')
    df.to_csv(path, index=False)
    card = DataCard(address='', file_type='csv', data_type='pandas.DataFrame')
    dataset_address = disk_datastore.upload_data('temp_dataframe.csv', path,
                                                 card)

    file_size = os.path.getsize(path)
    assert disk_datastore.get_file_size(dataset_address) == file_size


def test_disk_move_data(disk_datastore, alternate_disk_datastore, tmp_csv):
    # Put object
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')
    dataset_address = disk_datastore.upload_data('tmpcsv.csv', tmp_csv,
                                                 data_card)
    file_size = os.path.getsize(tmp_csv)
    assert disk_datastore.get_file_size(dataset_address) == file_size

    # Move object
    dest_address = 'deepchem://alternate-test/alternate-user/tmpcsv.csv'
    disk_datastore.move_object(dataset_address, dest_address,
                               alternate_disk_datastore)

    # Check that object is in new location
    assert disk_datastore.get_file_size(dataset_address) == 0
    assert alternate_disk_datastore.get_file_size(dest_address) == file_size


def test_disk_move_nested_data(disk_datastore, alternate_disk_datastore,
                               tmp_csv):
    # Put object
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')
    dataset_address = disk_datastore.upload_data('test move source/tmpcsv.csv',
                                                 tmp_csv, data_card)
    file_size = os.path.getsize(tmp_csv)
    assert disk_datastore.get_file_size(dataset_address) == file_size

    # Move object
    dest_address = 'deepchem://alternate-test/alternate-user/test move dest/tmpcsv.csv'
    disk_datastore.move_object(dataset_address, dest_address,
                               alternate_disk_datastore)

    # Check that object is in new location
    assert disk_datastore.get_file_size(dataset_address) == 0
    assert alternate_disk_datastore.get_file_size(dest_address) == file_size


def test_disk_copy_data(disk_datastore, alternate_disk_datastore, tmp_csv):
    # Put object
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')
    dataset_address = disk_datastore.upload_data('tmpcsv.csv', tmp_csv,
                                                 data_card)
    file_size = os.path.getsize(tmp_csv)
    assert disk_datastore.get_file_size(dataset_address) == file_size

    # Copy object
    dest_address = 'deepchem://alternate-test/alternate-user/tmpcsv_copy.csv'
    disk_datastore.copy_object(dataset_address, dest_address,
                               alternate_disk_datastore)

    # Check that object is copied to the new location
    assert alternate_disk_datastore.get_file_size(dest_address) == file_size


def test_disk_copy_nested_data(disk_datastore, alternate_disk_datastore,
                               tmp_csv):
    # Put object
    data_card = cards.DataCard(address='',
                               file_type='csv',
                               data_type='pandas.DataFrame',
                               description='this is a pandas dataframe')
    dataset_address = disk_datastore.upload_data('test copy source/tmpcsv.csv',
                                                 tmp_csv, data_card)
    file_size = os.path.getsize(tmp_csv)
    assert disk_datastore.get_file_size(dataset_address) == file_size

    # Copy object
    dest_address = 'deepchem://alternate-test/alternate-user/test copy dest/tmpcsv_copy.csv'
    disk_datastore.copy_object(dataset_address, dest_address,
                               alternate_disk_datastore)

    # Check that object is copied to the new location
    assert alternate_disk_datastore.get_file_size(dest_address) == file_size


def test_dir_disk_add_list_data(disk_datastore):
    """Test adding a directory to disk and list data"""
    disk_datastore.add_dir("test_dir/")
    assert "deepchem://" + disk_datastore.storage_loc + "/test_dir/" in disk_datastore.list_data(
    ).split("\n")


def test_dir_disk_delete(disk_datastore):
    """Test deleting a directory from disk."""
    disk_datastore.add_dir("test_dir_del/")
    assert "test/user/test_dir_del/" in disk_datastore.list_data()

    disk_datastore.delete_object("deepchem://test/user/test_dir_del/", 'dir')
    assert "test/user/test_dir_del/" not in disk_datastore.list_data()


def test_dir_disk_move(disk_datastore, alternate_disk_datastore):
    """Test moving a directory on disk."""
    disk_datastore.add_dir("test_dir_move/")
    assert "test/user/test_dir_move/" in disk_datastore.list_data()

    disk_datastore.move_object(
        "deepchem://test/user/test_dir_move/",
        "deepchem://alternate-test/alternate-user/test_dir_move/",
        alternate_disk_datastore)
    assert "test/user/test_dir_move/" not in disk_datastore.list_data()
    assert "alternate-test/alternate-user/test_dir_move/" in alternate_disk_datastore.list_data(
    )


def test_dir_disk_copy(disk_datastore, alternate_disk_datastore):
    """Test copying a directory on disk."""
    disk_datastore.add_dir("test_dir_copy/")
    assert "test/user/test_dir_copy/" in disk_datastore.list_data()

    disk_datastore.copy_object(
        "deepchem://test/user/test_dir_copy/",
        "deepchem://alternate-test/alternate-user/test_dir_copy/",
        alternate_disk_datastore)
    assert "alternate-test/alternate-user/test_dir_copy/" in alternate_disk_datastore.list_data(
    )


def test_disk_datastore_upload_model(disk_datastore):
    """Test uploading model."""
    model_card = cards.ModelCard(address='',
                                 model_type='random_forest_regressor',
                                 train_dataset_address='',
                                 description='this is a model')
    from sklearn.ensemble import RandomForestRegressor
    model = dc.models.SklearnModel(RandomForestRegressor())
    model.save()
    data_address = disk_datastore.upload_model('model', model, model_card)
    model_retrieved = disk_datastore.get_model(data_address)
    assert isinstance(model_retrieved, dc.models.Model)


def test_disk_datastore_upload_nested_model(disk_datastore):
    """Test uploading a nested model."""
    model_card = cards.ModelCard(address='',
                                 model_type='random_forest_regressor',
                                 train_dataset_address='',
                                 description='this is a model')
    from sklearn.ensemble import RandomForestRegressor
    model = dc.models.SklearnModel(RandomForestRegressor())
    model.save()
    data_address = disk_datastore.upload_model('test upload model/model', model,
                                               model_card)
    model_retrieved = disk_datastore.get_model(data_address)
    assert isinstance(model_retrieved, dc.models.Model)
    assert data_address == "deepchem://test/user/test upload model/model"


def test_disk_datastore_upload_model_from_memory(disk_datastore):
    """Test uploading model."""
    filenames = ['checkpoint1.pt', 'config.yaml']
    files = [
        open(
            os.path.join(os.path.dirname(__file__),
                         './assets/models/gcn/checkpoint1.pt'), 'rb'),
        open(
            os.path.join(os.path.dirname(__file__),
                         './assets/models/gcn/config.yaml'), 'rb')
    ]
    card = ModelCard(address='', train_dataset_address='', model_type='gcn')
    model_address = disk_datastore.upload_model_from_memory(
        model_name='model',
        model_files=files,
        model_filenames=filenames,
        card=card)

    # check that uploaded model exists
    card_retrieved = disk_datastore.get_card(model_address, 'model')
    assert card_retrieved.model_type == 'gcn'
    assert card_retrieved.address == 'deepchem://test/user/model'

    assert model_address == 'deepchem://test/user/model'


def test_disk_datastore_upload_nested_model_from_memory(disk_datastore):
    """Test uploading nested model from memory."""
    filenames = ['checkpoint1.pt', 'config.yaml']
    files = [
        open(
            os.path.join(os.path.dirname(__file__),
                         './assets/models/gcn/checkpoint1.pt'), 'rb'),
        open(
            os.path.join(os.path.dirname(__file__),
                         './assets/models/gcn/config.yaml'), 'rb')
    ]
    card = ModelCard(address='', train_dataset_address='', model_type='gcn')
    model_address = disk_datastore.upload_model_from_memory(
        model_name='test upload model from memory/model',
        model_files=files,
        model_filenames=filenames,
        card=card)

    # check that uploaded model exists
    card_retrieved = disk_datastore.get_card(model_address, 'model')
    assert card_retrieved.model_type == 'gcn'
    assert card_retrieved.address == 'deepchem://test/user/test upload model from memory/model'

    assert model_address == 'deepchem://test/user/test upload model from memory/model'
