import argparse
import functools
import logging
import pandas as pd
import numpy as np
from os.path import join
from utils._tools import read_data

import torch
from torch import nn
from torch.distributions import constraints
import pyro.poutine as poutine

import pyro
import pyro.distributions as dist
from pyro.infer import SVI, JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.optim import ClippedAdam

logging.basicConfig(format='%(relativeCreated) 9d %(message)s', level=logging.INFO)
device = torch.device("cpu") 
# This is a fully generative model of a batch of documents.
# data is a [num_words_per_doc, num_documents] shaped array of word ids
# (specifically it is not a histogram). We assume in this simple example
# that all documents have the same number of words.
def model(data=None, args=None, batch_size=123):
# def model(data=None, args=None, batch_size=None):

    # Globals.
    with pyro.plate("topics", args.num_topics, dim=-1):
        # HACK: topic_weights (num_topics=8,)
        topic_weights = pyro.sample(
            "topic_weights", dist.Gamma(1. / args.num_topics, 1.))
        assert topic_weights.shape==(args.num_topics,)

        # NOTE: beta prior for topic words
        with pyro.plate("prior_words", args.num_words_per_doc):
            topic_words = pyro.sample("topic_words", dist.Beta(
                torch.ones(1)*0.01, torch.ones(1)*0.01))
            assert topic_words.shape==(args.num_words_per_doc, args.num_topics)  

    # Locals.
    with pyro.plate("documents", args.num_docs) as ind:
        if data is not None:
            with pyro.util.ignore_jit_warnings():
                assert data.shape == (args.num_words_per_doc, args.num_docs)
            data = data[:, ind]
        
        doc_topics = pyro.sample("doc_topics", dist.Dirichlet(topic_weights))
        assert doc_topics.shape==(ind.size(0),args.num_topics)
        with pyro.plate("words", args.num_words_per_doc):
            # The word_topics variable is marginalized out during inference,
            # achieved by specifying infer={"enumerate": "parallel"} and using
            # TraceEnum_ELBO for inference. Thus we can ignore this variable in
            # the guide.
            word_topics = pyro.sample("word_topics", dist.Categorical(doc_topics),
                                      infer={"enumerate": "parallel"})
            ## word_topics (num_words=64, num_docs=1000)
            # NOTE: bernoulli likelihood
            word_indexes=torch.arange(0,args.num_words_per_doc).unsqueeze(1).repeat(1, ind.size(0))
            # assert word_indexes.shape==(args.num_words_per_doc, args.num_docs)
            data = pyro.sample("doc_words", dist.Bernoulli(topic_words[word_indexes, word_topics]),
                    obs=data)
            # assert data.shape==(args.num_words_per_doc, args.num_docs)

    return topic_weights, topic_words, data


# We will use amortized inference of the local topic variables, achieved by a
# multi-layer perceptron. We'll wrap the guide in an nn.Module.
def make_predictor(args):
    layer_sizes = ([args.num_words] +
                   [int(s) for s in args.layer_sizes.split('-')] +
                   [args.num_topics])
    logging.info('Creating MLP with sizes {}'.format(layer_sizes))
    layers = []
    for in_size, out_size in zip(layer_sizes, layer_sizes[1:]):
        layer = nn.Linear(in_size, out_size)
        layer.weight.data.normal_(0, 0.001)
        layer.bias.data.normal_(0, 0.001)
        layers.append(layer)
        layers.append(nn.Sigmoid())
    layers.append(nn.Softmax(dim=-1))
    return nn.Sequential(*layers)

def parametrized_guide(predictor, data, args, batch_size=None, print_args=False):
    topic_weights_posterior = pyro.param(
            "topic_weights_posterior",
            lambda: torch.ones(args.num_topics),
            constraint=constraints.positive)
    assert topic_weights_posterior.shape==(args.num_topics,)

    # print("shape of topic words")
    topic_words_posterior_a0 = pyro.param(
            "topic_words_posterior_a0",
            lambda: torch.ones(args.num_words_per_doc, args.num_topics)*0.01,
            # lambda: torch.ones(args.num_topics)*0.5,
            constraint=constraints.positive)
    # assert topic_words_posterior_a0==（）
    topic_words_posterior_a1 = pyro.param(
            "topic_words_posterior_a1",
            lambda: torch.ones(args.num_words_per_doc, args.num_topics)*0.01,
            # lambda: torch.ones(args.num_topics)*0.5,
            constraint=constraints.positive)

    with pyro.plate("topics", args.num_topics):
        pyro.sample("topic_weights", dist.Gamma(topic_weights_posterior, 1.))
        # NOTE: phi distribution
        # pyro.sample("topic_words", dist.Dirichlet(topic_words_posterior))
        # with pyro.plate("words", args.num_words_per_doc):
        with pyro.plate("prior_words", args.num_words_per_doc):
            topic_words=pyro.sample("topic_words", dist.Beta(
                # torch.ones([1],dtype=torch.int64)*0.5, torch.tensor([1],dtype=torch.int64)*0.5))
                topic_words_posterior_a0, topic_words_posterior_a1))
            # assert topic_words.shape==(args.num_words_per_doc, args.num_topics)
        # pyro.sample("topic_words", dist.Beta(topic_words_posterior_a0,topic_words_posterior_a1))

    # Use an amortized guide for local variables.
    pyro.module("predictor", predictor)
    with pyro.plate("documents", args.num_docs, batch_size) as ind:
        data = data[:, ind].type(torch.int64)
        # The neural network will operate on histograms rather than word
        # index vectors, so we'll convert the raw data to a histogram.
        # print(data)
        # print(data.shape)
        # print(data.dtype)
        # quit()
        counts = (torch.zeros(args.num_words_per_doc, ind.size(0))
                       .scatter_add(0, data, torch.ones(data.shape)))
        doc_topics = predictor(counts.transpose(0, 1))
        pyro.sample("doc_topics", dist.Delta(doc_topics, event_dim=1))

    if(print_args):
        print("shape of final θ: ")
        print(doc_topics.shape)
        print(doc_topics)
        # print("sum of row:")
        # print(doc_topics.sum(axis=1))
        # print("sum of column: ")
        # print(doc_topics.sum(axis=0))
        print("shape of final φ: ")
        print(topic_words)
        print(topic_words.shape)
        # print("sum of row:")
        # print(topic_words.sum(axis=1))
        # print("sum of column: ")
        # print(topic_words.sum(axis=0))


def main(args):
    logging.info('Generating data')
    pyro.set_rng_seed(0)
    pyro.clear_param_store()
    pyro.enable_validation(__debug__)

    rawdata =torch.tensor(read_data(
        "/data/liu/mimic3/CLAMP_NER/single_drug_analysis/FEATURE/PRE_PROCESS/pres_rxnorm_matrix.csv").\
            set_index("HADM_ID").values.T, dtype=torch.int64)
    # print(list(filter(lambda x: x<2, data.sum(axis=0))))
    # rawdata=None
    # data.sum(axis=0)
    # quit()
    true_topic_weights, true_topic_words, data = model(args=args,data=rawdata)
    print(data.shape)
    # quit()

    # quit()
    # print(data.shape)
    # quit()
    # data = torch.utils.data.DataLoader(data_raw,batch_size=args.batch_size,shuffle=False)
    # We can generate synthetic data directly by calling the model.
    # true_topic_weights, true_topic_words, data = model(args=args)
    # print(data)
    # print(data.shape)
    # print("\ntopic_weights_prior:")
    # print(true_topic_weights.shape)
    # print("\nshape of prior φ:")
    # print(true_topic_words.shape)
    # quit()

    # We'll train using SVI.
    logging.info('-' * 40)
    logging.info('Training on {} documents'.format(args.num_docs))
    predictor = make_predictor(args)
    guide = functools.partial(parametrized_guide, predictor)
    Elbo = JitTraceEnum_ELBO if args.jit else TraceEnum_ELBO
    elbo = Elbo(max_plate_nesting=2)
    optim = ClippedAdam({'lr': args.learning_rate})
    svi = SVI(model, guide, optim, elbo)
    logging.info('Step\tLoss')
    for step in range(args.num_steps):
        loss = svi.step(data, args=args, batch_size=args.batch_size)
        if step % 10 == 0:
            logging.info('{: >5d}\t{}'.format(step, loss))
    guide = functools.partial(guide,  print_args=True)
    loss = elbo.loss(model, guide, data, args=args)
    logging.info('final loss = {}'.format(loss))


if __name__ == '__main__':
    

    assert pyro.__version__.startswith('1.5.2')
    parser = argparse.ArgumentParser(description="Amortized Latent Dirichlet Allocation")
    parser.add_argument("-t", "--num-topics", default=8, type=int)
    parser.add_argument("-w", "--num-words", default=1624, type=int)
    parser.add_argument("-d", "--num-docs", default=49955, type=int)
    # parser.add_argument("-d", "--num-docs", default=32, type=int)
    # parser.add_argument("-wd", "--num-words-per-doc", default=12, type=int)
    parser.add_argument("-wd", "--num-words-per-doc", default=1624, type=int)
    parser.add_argument("-n", "--num-steps", default=1000, type=int)
    parser.add_argument("-l", "--layer-sizes", default="100-100")
    parser.add_argument("-lr", "--learning-rate", default=0.01, type=float)
    parser.add_argument("-b", "--batch-size", default=32, type=int)
    parser.add_argument('--jit', action='store_true')
    args = parser.parse_args()
    main(args)



