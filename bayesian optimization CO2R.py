import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
import pickle


#%% plot activity & selectivity map
with open('.../activity map and selectivity map.pkl', 'rb') as f: 
    x_map_activity, y_map_activity, z_map_activity, x_map_selectivity, y_map_selectivity, z_map_selectivity = pickle.load(f)


interp_activity = interpolate.RegularGridInterpolator((x_map_activity, y_map_activity), z_map_activity,
                                  bounds_error=False, fill_value=None)

xnew = np.arange(-0.6, 0.8, 0.01)
ynew = np.arange(-2, 0.5, 0.01)
xxnew, yynew = np.meshgrid(xnew, ynew, indexing='xy')
znew = interp_activity((xxnew, yynew))

fig, axes = plt.subplots(figsize=(4, 3))
im = plt.imshow(znew, extent=[min(xnew),max(xnew),min(ynew),max(ynew)], origin="lower", cmap=cm.RdBu_r, alpha=0.5, aspect=.4)
plt.colorbar(im,fraction=0.034, pad=0.04)
axes.set_xlabel('H adsorption energy (eV)', fontsize = 12, fontname="Arial")
axes.set_ylabel('CO adsorption energy (eV)', fontsize = 12, fontname="Arial")
plt.show()


interp_selectivity = interpolate.RegularGridInterpolator((x_map_selectivity, y_map_selectivity), z_map_selectivity,
                                  bounds_error=False, fill_value=None)
xnew = np.arange(-0.6, 0.8, 0.01)
ynew = np.arange(-2, 0.5, 0.01)
znew = interp_selectivity((xxnew, yynew))

fig, axes = plt.subplots(figsize=(4, 3))
im = plt.imshow(znew, extent=[min(xnew),max(xnew),min(ynew),max(ynew)], origin="lower", cmap=cm.RdBu_r, alpha=0.5, aspect=.4)
plt.colorbar(im,fraction=0.034, pad=0.04)
axes.set_xlabel('H adsorption energy (eV)', fontsize = 12, fontname="Arial")
axes.set_ylabel('CO adsorption energy (eV)', fontsize = 12, fontname="Arial")
plt.show()

#%% caluclate activity & selectivity

with open('.../dataset.pkl', 'rb') as f:
    df = pickle.load(f)


deltaE_H = df.deltaE_H.to_numpy()
deltaE_CO = df.deltaE_CO.to_numpy()

a = np.logical_and(deltaE_H>min(x_map_activity), deltaE_H<max(x_map_activity))
b = np.logical_and(deltaE_CO>min(y_map_activity), deltaE_CO<max(y_map_activity))
idx = np.logical_and(a, b)

activity = np.empty(len(deltaE_H))
activity[:] = np.nan
selectivity = np.empty(len(deltaE_H))
selectivity[:] = np.nan

for i in np.where(idx==True)[0]:
    activity[i] = interpolate.interpn((x_map_activity, y_map_activity), z_map_activity, np.array([deltaE_H[i], deltaE_CO[i]]))
    selectivity[i] = interpolate.interpn((x_map_selectivity, y_map_selectivity), z_map_selectivity, np.array([deltaE_H[i], deltaE_CO[i]]))


df['activity'] = activity
df['selectivity'] = selectivity

df = df.drop(np.where(idx==False)[0])
df.index = range(len(df))


#%% inputs
chemicals_unique = []
for i in range(len(df)):
    chemicals_unique.append(list(set(df.iloc[i].POSCAR_CO.get_chemical_symbols())))

flat_chemicals_unique = [item for sublist in chemicals_unique for item in sublist]
chemicals_all_unique = list(set(flat_chemicals_unique))

# inputs - POSCAR_H
chemicals_list = []
positions_list = []
for i in range(len(df)):
    chemicals_list.append(df.iloc[i].POSCAR_H.get_chemical_symbols())
    positions_list.append(df.iloc[i].POSCAR_H.get_positions().tolist())

# normalized position
positions_min = min(min(min(positions_list)))
positions_max = max(max(max(positions_list)))
positions_list_norm = []
for position in positions_list:
    position_array = np.asarray(position)
    position_array = (position_array- positions_min) / (positions_max - positions_min)
    position = position_array.tolist()
    positions_list_norm.append(position)

# atom one hot
chemicals_onehot_list = []
for chemical in chemicals_list:
    chemical_onehot = []
    for atom in chemical:
        chemical_onehot.append([int(element == atom) for element in chemicals_all_unique])
    chemicals_onehot_list.append(chemical_onehot)

# atom descriptor
periodic_table = pd.read_excel('.../List of elements by atomic properties.xlsx')

chemicals_descriptor_list = []
for chemical in chemicals_list:
    chemical_descriptor = []
    for atom in chemical:
        chemical_descriptor.append(periodic_table.loc[periodic_table['symbol']==atom, ['atomic mass norm', 'electronegativity norm', 'atomic radii norm']].values[0].tolist()) # ['atomic mass', 'electronegativity', 'atomic radii'] 'number', 
    chemicals_descriptor_list.append(chemical_descriptor)

poscar_list = [np.hstack((x,y,z)).tolist() for x,y,z in zip(chemicals_onehot_list, chemicals_descriptor_list, positions_list_norm)]

for sublist in poscar_list:
    sublist[:] = sublist + [np.zeros(len(chemicals_onehot_list[0][0])+len(chemicals_descriptor_list[0][0])).tolist()+[0.,0.,0.]] * (73 - len(sublist))

poscar_N = np.asarray(poscar_list)


# inputs - POSCAR_CO
chemicals_list = []
positions_list = []
for i in range(len(df)):
    chemicals_list.append(df.iloc[i].POSCAR_CO.get_chemical_symbols())
    positions_list.append(df.iloc[i].POSCAR_CO.get_positions().tolist())


# normalized position
positions_min = min(min(min(positions_list)))
positions_max = max(max(max(positions_list)))
positions_list_norm = []
for position in positions_list:
    position_array = np.asarray(position)
    position_array = (position_array- positions_min) / (positions_max - positions_min)
    position = position_array.tolist()
    positions_list_norm.append(position)


# atom one hot
chemicals_onehot_list = []
for chemical in chemicals_list:
    chemical_onehot = []
    for atom in chemical:
        chemical_onehot.append([int(element == atom) for element in chemicals_all_unique])
    chemicals_onehot_list.append(chemical_onehot)


# atom descriptor
chemicals_descriptor_list = []
for chemical in chemicals_list:
    chemical_descriptor = []
    for atom in chemical:
        chemical_descriptor.append(periodic_table.loc[periodic_table['symbol']==atom, ['atomic mass norm', 'electronegativity norm', 'atomic radii norm']].values[0].tolist()) #['atomic mass', 'electronegativity', 'atomic radii'] 'number', 
    chemicals_descriptor_list.append(chemical_descriptor)

poscar_list = [np.hstack((x,y,z)).tolist() for x,y,z in zip(chemicals_onehot_list, chemicals_descriptor_list, positions_list_norm)]

for sublist in poscar_list:
    sublist[:] = sublist + [np.zeros(len(chemicals_onehot_list[0][0])+len(chemicals_descriptor_list[0][0])).tolist()+[0.,0.,0.]] * (74 - len(sublist))


poscar_NOH = np.asarray(poscar_list)

# inputs
poscar = np.concatenate((poscar_N, poscar_NOH), axis=1)


#%% The deterministic model
import tensorflow as tf

# Define model
class DeepResNet(tf.keras.Model):
  """Defines a multi-layer residual network."""
  def __init__(self, num_classes, kernel_size=2, num_conv2d_layers=3, num_dense_layers=3,
               num_conv2d_hidden=64, num_dense_hidden=64, dropout_rate=0.1, **classifier_kwargs):
    super().__init__()
    # Defines class meta data.
    self.kernel_size = kernel_size
    self.num_conv2d_hidden = num_conv2d_hidden
    self.num_dense_hidden = num_dense_hidden
    self.num_conv2d_layers = num_conv2d_layers
    self.num_dense_layers = num_dense_layers
    self.dropout_rate = dropout_rate
    self.classifier_kwargs = classifier_kwargs
    
    # Defines the hidden layers.
    self.input_layer = tf.keras.layers.Conv2D(filters=self.num_conv2d_hidden, kernel_size=1, trainable=False)
    self.conv2d_layers = [self.make_conv2d_layer() for _ in range(num_conv2d_layers)]
    
    # Defines the output layer.
    self.classifier = self.make_output_layer(num_classes)
    
  def call(self, inputs):
    # Projects the 2d input data to high dimension.
    hidden = self.input_layer(inputs)
    
    # Computes the ResNet hidden representations.
    for i in range(self.num_conv2d_layers):
      resid = self.conv2d_layers[i](hidden)
      hidden += resid
    
    hidden = tf.keras.layers.GlobalMaxPooling2D()(hidden)
    hidden = tf.keras.layers.Flatten()(hidden)
      
    return self.classifier(hidden)
    
  def make_conv2d_layer(self):
    """Uses the Conv2d layer as the hidden layer."""
    return tf.keras.layers.Conv2D(filters=self.num_conv2d_hidden, kernel_size=self.kernel_size, strides=1, padding='same', activation='relu')
    
  def make_output_layer(self, num_classes):
    """Uses the Dense layer as the output layer."""
    return tf.keras.layers.Dense(
        num_classes, **self.classifier_kwargs)


# %% SNGP

import sys
sys.path.append(".../layers")
import gaussian_process
import spectral_normalization


class DeepResNetSNGP(DeepResNet):
  def __init__(self, spec_norm_bound=0.95, **kwargs):
    self.spec_norm_bound = spec_norm_bound
    super().__init__(**kwargs)
    
  def make_conv2d_layer(self):
    """Applies spectral normalization to the conv2d hidden layer."""
    conv2d_layer = super().make_conv2d_layer()
    return spectral_normalization.SpectralNormalizationConv2D( 
        conv2d_layer, norm_multiplier=self.spec_norm_bound)
    
  def make_dense_layer(self):
    """Applies spectral normalization to the dense hidden layer."""
    dense_layer = super().make_dense_layer()
    return spectral_normalization.SpectralNormalization( 
        dense_layer, norm_multiplier=self.spec_norm_bound)
    
  def make_output_layer(self, num_classes):
    """Uses Gaussian process as the output layer."""
    return gaussian_process.RandomFeatureGaussianProcess( 
        num_classes,
        gp_cov_momentum=-1,
        **self.classifier_kwargs)
    
  def call(self, inputs, training=False, return_covmat=False):
    # Gets logits and a covariance matrix from the GP layer.
    logits, covmat = super().call(inputs)
    
    # Returns only logits during training.
    if not training and return_covmat:
      return logits, covmat
      
    return logits



class ResetCovarianceCallback(tf.keras.callbacks.Callback):
  
  def on_epoch_begin(self, epoch, logs=None):
    """Resets covariance matrix at the beginning of the epoch."""
    if epoch > 0:
      self.model.classifier.reset_covariance_matrix()

class DeepResNetSNGPWithCovReset(DeepResNetSNGP):
  def fit(self, *args, **kwargs):
    """Adds ResetCovarianceCallback to model callbacks."""
    kwargs["callbacks"] = list(kwargs.get("callbacks", []))
    kwargs["callbacks"].append(ResetCovarianceCallback())
    
    return super().fit(*args, **kwargs)

#%% Constrained Bayesian Optimization
    
# select initial samples 
import random
symbol = df['comp'].to_numpy()
remaining_idx = list(range(len(df)))
training_idx = []

symbol_unique, symbol_counts = np.unique(symbol, return_counts=True)
symbol_size = dict(zip(symbol_unique, symbol_counts))
symbol_idx = []
for i in range(len(symbol_unique)):
    symbol_idx.append(np.where(symbol==symbol_unique[i])[0].tolist())

for i in range(len(symbol_unique)):
    symbol_list = symbol_idx[i]
    if len(symbol_list) > 11:
        symbol_list_idx = random.sample(symbol_list, 1)
        training_idx.append(symbol_list_idx)

training_idx = [item for sublist in training_idx for item in sublist]

remaining_idx = list(range(len(df)))
remaining_idx = list(set(remaining_idx) - set(training_idx))



#% record training data
training_idx_all = []
training_idx_all.append(training_idx)


constraint = (df.selectivity.to_numpy()>0.9) + 0
for iteration in range(80):
    # regression
    # data
    X_train = poscar[training_idx]
    y_train = df.activity.to_numpy()[training_idx]
    X_test = poscar[remaining_idx]
    y_test = df.activity.to_numpy()[remaining_idx]
    
    # train model
    resnet_config = dict(num_classes=1, kernel_size=(1, X_train.shape[2]),
                         num_conv2d_layers=5, num_dense_layers=2, 
                         num_conv2d_hidden=64, num_dense_hidden=128)
    
    loss = tf.keras.losses.MeanSquaredError()
    metrics = tf.keras.losses.MeanSquaredError()
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)
    train_config = dict(loss=loss, metrics=metrics, optimizer=optimizer)
    fit_config = dict(batch_size=32, epochs=500)

    sngp_model_regression = DeepResNetSNGPWithCovReset(**resnet_config)
    sngp_model_regression.compile(**train_config)
    
    sngp_model_regression.fit(X_train[:, np.newaxis, :, :], y_train, **fit_config)
    
    
    # prediction
    sngp_logits, sngp_covmat = sngp_model_regression(X_test[:, np.newaxis, :, :], return_covmat=True)
    sngp_variance = tf.linalg.diag_part(sngp_covmat)[:, None]
    y_pred = sngp_logits.numpy().flatten()
    y_pred_std = np.sqrt(sngp_variance.numpy()).flatten()
    
    
    # expected improvement
    from scipy.stats import norm
    if np.sum(constraint[training_idx]>0) > 0:
        ft_min = -np.max(y_train[constraint[training_idx]>0])
    else:
        ft_min = 0
    mu_x = -y_pred
    sigma_x = y_pred_std
    EI = (ft_min - mu_x) * norm.cdf((ft_min - mu_x)/sigma_x) + sigma_x * norm.pdf((ft_min - mu_x)/sigma_x)
    
    
    # classification
    # data
    X_train = poscar[training_idx]
    y_train = constraint[training_idx]
    X_test = poscar[remaining_idx]
    y_test = constraint[remaining_idx]
    
    neg = np.sum(y_train==0)
    pos = np.sum(y_train==1)
    total = len(X_train)
    weight_for_0 = (1 / neg) * (total / 2.0)
    weight_for_1 = (1 / pos) * (total / 2.0)
    class_weight = {0: weight_for_0, 1: weight_for_1}

    # train model
    resnet_config = dict(num_classes=2, kernel_size=(1, X_train.shape[2]),
                         num_conv2d_layers=5, num_dense_layers=2, 
                         num_conv2d_hidden=64, num_dense_hidden=128)
    
    loss = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
    metrics = tf.keras.metrics.SparseCategoricalAccuracy(),
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-5)
    train_config = dict(loss=loss, metrics=metrics, optimizer=optimizer)
    fit_config = dict(batch_size=32, epochs=300, class_weight=class_weight)

    sngp_model_classification = DeepResNetSNGPWithCovReset(**resnet_config)
    sngp_model_classification.compile(**train_config)
    
    sngp_model_classification.fit(X_train[:, np.newaxis, :, :], y_train, **fit_config)
    
    # prediction
    sngp_logits, sngp_covmat = sngp_model_classification(X_test[:, np.newaxis, :, :], return_covmat=True)
    sngp_variance = tf.linalg.diag_part(sngp_covmat)[:, None]
    sngp_logits_adjusted = sngp_logits / tf.sqrt(1. + (np.pi / 8.) * sngp_variance) # Eq. 19
    sngp_probs = tf.nn.softmax(sngp_logits_adjusted, axis=-1)[:, 1]
    uncertainty = sngp_probs * (1. - sngp_probs)
    uncertainty = (uncertainty - np.min(uncertainty)) / (np.max(uncertainty) - np.min(uncertainty))
    
    
    # constrained expected improvement
    PF = sngp_probs
    EIC = PF * EI
    
    # next samples
    k_next = 1
    EIC_idx = np.argpartition(EIC, -k_next)
    next_idx = np.array(remaining_idx)[EIC_idx[-k_next:]].tolist()
    
    df.activity.to_numpy()[next_idx]
    constraint[next_idx]
    
    training_idx_all.append(next_idx)
    training_idx = [item for sublist in training_idx_all for item in sublist]
    remaining_idx = list(set(remaining_idx) - set(next_idx))
    
# end of iteration #


#%% plot constrained BO sampling
from matplotlib.lines import Line2D

total_itr = 81

constraint = (df.selectivity.to_numpy() > 0.999) + 0

xnew = np.arange(-0.6, 0.8, 0.01)
ynew = np.arange(-2, 0.5, 0.01)
xxnew, yynew = np.meshgrid(xnew, ynew, indexing='xy')
znew = interp_activity((xxnew, yynew))

fig, axes = plt.subplots(figsize=(4, 3))
im = plt.imshow(znew, extent=[min(xnew),max(xnew),min(ynew),max(ynew)], origin="lower", cmap=cm.RdBu_r, alpha=0.5, aspect=.4)
plt.colorbar(im,fraction=0.034, pad=0.04)
axes.set_xlabel('H adsorption energy (eV)', fontsize = 12, fontname="Arial")
axes.set_ylabel('CO adsorption energy (eV)', fontsize = 12, fontname="Arial")


for iteration in range(total_itr):
    legend = []
    for itr in [iteration]: #range(0)
        fc = np.where(constraint[training_idx_all[itr]] == 0, 'none', 'C'+str(itr))
        ec = np.where(constraint[training_idx_all[itr]] == 0, 'C'+str(itr), 'white')
        pt = plt.scatter(df.deltaE_H.to_numpy()[training_idx_all[itr]],
                    df.deltaE_CO.to_numpy()[training_idx_all[itr]],
                    facecolors=fc,
                    edgecolors=ec,
                    linewidths=1)
        point0 = Line2D([0], [0], label='Itr %s' % itr + '_0', marker='o',
                 markeredgecolor='C'+str(itr), markerfacecolor='none', linestyle='')
        point1 = Line2D([0], [0], label='Itr %s' % itr + '_1', marker='o',
                 markeredgecolor='white', markerfacecolor='C'+str(itr), linestyle='')
        legend.append(point0)
        legend.append(point1)
plt.show()

# top 10
idx_bo = [x for y in training_idx_all[:total_itr] for x in y]
idx_selected = np.array(idx_bo)[constraint[idx_bo] == 1]
np.sort(df['activity'].to_numpy()[idx_selected])
print(df[['activity', 'selectivity', 'comp', 'mpid', 'miller', 'deltaE_H', 'deltaE_CO']].loc[idx_selected].sort_values(by='activity', ascending=False)[:10])