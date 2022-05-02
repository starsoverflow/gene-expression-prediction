import timeit

import numpy
import numpy as np
import torch
import collections
import pdb
import torch.utils.data
import csv
import json

from sklearn import metrics
from torch import nn
import torch.nn.functional as F
from torchvision import transforms, utils
import math

import data_loader
from data_loader import SingleCellData


def compute_aupr(all_targets, all_predictions):

    precision, recall, thresholds = metrics.precision_recall_curve(all_targets, all_predictions,
                                                                   pos_label=1)
    auPR = metrics.auc(recall, precision)  # ,reorder=True)

    return auPR


def run_one_epoch(train_flag, dataloader, model, optimizer, device="cpu"):
    torch.set_grad_enabled(train_flag)
    model.train() if train_flag else model.eval()

    losses = []
    accuracies = []
    auprs = []

    for (x,y) in dataloader:  # collection of tuples with iterator

        x = x.to(device)
        y = y.to(device).to(torch.float32)

        output = model(x)  # forward pass
        output = output.squeeze()  # remove spurious channel dimension

        loss = F.binary_cross_entropy_with_logits( output, y ) # numerically stable

        if train_flag:
            loss.backward() # back propagation
            optimizer.step()
            optimizer.zero_grad()

        losses.append(loss.detach().cpu().numpy())
        accuracy = torch.mean( ( (output > .5) == (y > .5) ).float() )
        accuracies.append(accuracy.detach().cpu().numpy())

        aupr = compute_aupr(y.cpu().numpy(), output.cpu().detach().numpy())
        auprs.append(aupr)

    return (np.mean(losses), np.mean(accuracies), np.mean(auprs))


class CNN_LSTM_1d(nn.Module):

    def __init__(self, seq_len, n_in_channels=1,
                 n_output_channels=1,
                 dropout=0.25):

        super(CNN_LSTM_1d, self).__init__()

        self.seq_len = seq_len  # amount of sequence context required
        conv_output_len = self.seq_len

        conv_layers = []

        for i in range(2):
            in_channels = 32
            out_channels = 32
            kernel_size = 3
            if i == 0:
                in_channels = n_in_channels
                kernel_size = 15
            conv_layers += [nn.Conv1d(in_channels=in_channels, out_channels=out_channels, kernel_size=kernel_size),
                            nn.BatchNorm1d(out_channels),
                            nn.ReLU(inplace=True),
                            nn.Dropout(dropout),
                            ]
            conv_output_len = (conv_output_len - kernel_size) + 1
            conv_output_len = conv_output_len
        print("Sequence length after convolution:", conv_output_len)

        # If you have a model with lots of layers, you can create a list first and
        # then use the * operator to expand the list into positional arguments, like this:
        self.conv_net = nn.Sequential(*conv_layers)

        rnn_hidden = 256
        self.rnn_net = nn.LSTM(32, rnn_hidden, 2, bidirectional=True, dropout=0.5)

        self.dense_net = nn.Sequential(nn.Dropout(0.5),
                                       nn.Linear(rnn_hidden * conv_output_len * 2, 256),
                                       nn.Dropout(dropout),
                                       nn.ELU(inplace=True),
                                       nn.Linear(256, 32),
                                       nn.Dropout(dropout),
                                       nn.ELU(inplace=True),
                                       nn.Linear(32, n_output_channels))

    def forward(self, x):
        net = self.conv_net(x)

        # print(net.shape)
        net = net.transpose(1, 2)
        net = net.transpose(0, 1).contiguous()
        # print(net.shape)

        lstm_out, _ = self.rnn_net(net)
        # print(lstm_out.shape)

        net = lstm_out.transpose(0, 1).contiguous()
        # print(net.shape)
        net = net.view(net.size(0), -1)
        # print(net.shape)
        net = self.dense_net(net)
        return (net)


def train_model(model, train_dataloader, validation_dataloader, epochs=100, patience=10, verbose=True):
    """
    Train a 1D CNN model and record accuracy metrics.
    """
    device = "cpu"
    model.to(device)

    optimizer = torch.optim.Adam(model.parameters(), amsgrad=True, weight_decay=0.001)

    # 3. Run the training loop with early stopping.
    train_accs = []
    val_accs = []
    patience_counter = patience
    best_val_loss = np.inf
    check_point_filename = 'model_checkpoint.pt'  # to save the best model fit to date
    for epoch in range(epochs):
        start_time = timeit.default_timer()
        train_loss, train_acc, train_aupr = run_one_epoch(True, train_dataloader, model, optimizer, device)
        val_loss, val_acc, val_aupr = run_one_epoch(False, validation_dataloader, model, optimizer, device)
        train_accs.append(train_acc)
        val_accs.append(val_acc)
        if val_loss < best_val_loss:
            torch.save(model.state_dict(), check_point_filename)
            best_val_loss = val_loss
            patience_counter = patience
        else:
            patience_counter -= 1
            if patience_counter <= 0:
                model.load_state_dict(torch.load(check_point_filename))  # recover the best model so far
                break
        elapsed = float(timeit.default_timer() - start_time)
        if verbose:
            print("Epoch %i took %.2fs. Train loss: %.4f acc: %.4f aupr %.4f. Val loss: %.4f acc: %.4f aupr %.4f. Patience left: %i" %
                  (epoch + 1, elapsed, train_loss, train_acc, train_aupr, val_loss, val_acc, val_aupr, patience_counter))

    _, final_train_acc, final_train_aupr = run_one_epoch(False, train_dataloader, model, optimizer, device)
    _, final_val_acc, final_val_aupr = run_one_epoch(False, validation_dataloader, model, optimizer, device)
    _, final_test_acc, final_test_aupr = run_one_epoch(False, test_dataloader, model, optimizer, device)
    print("Train acc: %.4f aupr %.4f" % (final_test_acc, final_test_aupr))
    return model, train_accs, val_accs, final_train_acc, final_val_acc


train_dataloader, val_dataloader, test_dataloader, channels, context_length = data_loader.load_data()

my_cnn_lstm_1d = CNN_LSTM_1d(context_length, 2, 1)

print("Input length:", my_cnn_lstm_1d.seq_len)


my_cnn_lstm_1d, train_accs, val_accs, _, _ = train_model(my_cnn_lstm_1d, train_dataloader, val_dataloader,
                                                         epochs=100, patience=20)


