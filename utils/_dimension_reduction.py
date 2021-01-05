import tensorflow as tf
import keras

from keras.backend.tensorflow_backend import set_session
from keras.models import Sequential, Model
from keras import layers
from keras import optimizers
from keras.layers import  BatchNormalization, Activation,  Dense, Input, Dropout
# from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import MinMaxScaler

from os.path import join
from _path import *
from _tools import *

# NOTE: Test whether tensorflow can work
# print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))


def reduce_dimension(
    inputdata, index_col, file_name="dissum_Autoencoder", 
    new_dimension=128, args=["sigmoid","binary_crossentropy"],
    data_transformed=False):
    physical_devices = tf.config.list_physical_devices('GPU') 
    for gpu_instance in physical_devices: 
        tf.config.experimental.set_memory_growth(gpu_instance, True)  

    num_dimension=inputdata.shape[1]

    ## Build autoencoder model
    model = Sequential()
    model.add(Dense(int(num_dimension/2),  activation='relu', input_shape=(num_dimension,)))
    model.add(BatchNormalization())
    # model.add(Dense(800, activation='relu'))
    model.add(Dense(new_dimension,  activation='relu', name="bottleneck"))
    model.add(BatchNormalization())
    model.add(Dense(int(num_dimension/2), activation='relu'))
    # model.add(Dense(1600,  activation='relu'))
    model.add(Dense(num_dimension,  activation=args[0]))
    model.summary()

    encoder = Model(model.input, model.get_layer('bottleneck').output)
    model.compile(loss=args[1], optimizer = optimizers.Adadelta())
    # model.compile(loss='binary_crossentropy', optimizer=optimizers.adam)
    if(data_transformed):
        scaler = MinMaxScaler()
        inputdata= scaler.fit_transform(inputdata)

    model.fit(inputdata,inputdata, batch_size=new_dimension, epochs=500,shuffle=True)
    repre_128 = encoder.predict(inputdata)
    repre_128_df = pd.DataFrame(repre_128)
    # print(repre_128_df.head())
    write2file(repre_128_df,join(singledrug_feature_prefix,"%s_128"%file_name))
    write2file(pd.DataFrame({"HADM_ID":index_col}),join(singledrug_feature_prefix,"%s_EPIS"%file_name))



# def preprocess_inputdata( inputdatum):
#     for i in inputdatum:
#         print(i)

def get_dissum_repre128():
    file_list=[
        join(singledrug_featurepreprocess_prefix,"diag_matrix"),
        join(singledrug_featurepreprocess_prefix,"pres_matrix"),
        join(concat_clamp_prefix,"allepis_newICD9_CODE"),
        join(concat_clamp_prefix,"allepis_newRxNorm")]
    
    data=[read_data(file,dtype={"HADM_ID":str}) for file in file_list]
    [print(da.shape) for da in data]
    data=[df.set_index("HADM_ID") for df in data]
    all_data=pd.concat(data,axis=1,sort=True).fillna(0)
    print(all_data.shape)
    print(all_data.loc[:, (all_data!= 0).any(axis=0)].shape)

    # reduce_dimension(all_data.loc[:, (all_data!= 0).any(axis=0)],list(all_data.index.values))


def get_mimic_five_repre128():
    file_list=[
        join(singledrug_featurepreprocess_prefix,"diag_matrix"),
        join(singledrug_featurepreprocess_prefix,"pres_matrix"),
        join(singledrug_featurepreprocess_prefix,"procedure_matrix"),
        join(singledrug_featurepreprocess_prefix,"lab_matrix")]

    
    data=[read_data(file,dtype={"HADM_ID":str}) for file in file_list]
    [print(datadf.shape) for datadf in data]
    # all_data=data[0]
    # print(all_data['HADM_ID'])

    # for data_ in data[1:]:
    #     all_data = inner_join(all_data,data_,"HADM_ID")
    
    data=[df.set_index("HADM_ID") for df in data]
    all_data=pd.concat(data,axis=1,sort=True).fillna(0)
    print(all_data.shape)
    all_data = all_data.rename_axis('HADM_ID').reset_index()


    demo_df = read_data(join(singledrug_featurepreprocess_prefix,"demographic_matrix"))
    patient_epis_df = read_data(join(singledrug_featurepreprocess_prefix,"patient_epis_map"),dtype={"HADM_ID":str})
    epis_demo_df = left_join(patient_epis_df,demo_df,"SUBJECT_ID")
    print(demo_df.shape)
    

    all_data=inner_join(all_data,epis_demo_df,"HADM_ID").drop('SUBJECT_ID', 1).fillna(0).set_index("HADM_ID")
    print(all_data.shape)
    all_data=all_data.loc[:, (all_data!= 0).any(axis=0)]
    print(all_data.shape)

    reduce_dimension(
        all_data,
        list(all_data.index.values),
        file_name="five_Autoencoder",
        # args=["sigmoid","mean_squared_error"],
        data_transformed=True)

# get_mimic_five_repre128()
# get_dissum_repre128()

# preprocess_inputdata("aa", "bb", "dd", "ee")


