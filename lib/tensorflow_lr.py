"""Tensorflow implementation of a binary logistic classification."""
import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt

def variable_summaries(var,name):
  """Attach a lot of summaries to a Tensor (for TensorBoard visualization)."""
  with tf.name_scope(name):
    mean = tf.reduce_mean(var)
    tf.summary.scalar('mean', mean)
    with tf.name_scope('stddev'):
      stddev = tf.sqrt(tf.reduce_mean(tf.square(var - mean)))
    tf.summary.scalar('stddev', stddev)
    tf.summary.scalar('max', tf.reduce_max(var))
    tf.summary.scalar('min', tf.reduce_min(var))
    tf.summary.histogram('histogram', var)


class LR():
    def __init__(self,train_set,val_set,test_set,learning_rate=0.0005,epochs=100):
        self.alpha, self.epochs = learning_rate, epochs

        self.X_train = train_set[0]
        self.y_train = train_set[1]

        self.X_val = val_set[0]
        self.y_val = val_set[1]

        self.X_test = test_set[0]
        self.y_test = test_set[1]

        # initialize n equal to the number of columns in the training dataframe
        n = self.X_train.shape[1]

        # placeholder with length equal to the number of features
        self.X = tf.placeholder(tf.float32,[None,n])

        # placeholder for Y which can only hold 2 values since it is a 2 class classification
        self.Y = tf.placeholder(tf.float32,[None,2])

        # trainable variable weights
        self.W = tf.Variable(tf.zeros([n,2]))
        variable_summaries(self.W,'Weights')

        # trainable variable bias
        self.b = tf.Variable(tf.zeros([2]))
        variable_summaries(self.b,'Bias')

        # hypothesis for logistic regression with a sigmoid activation function
        self.Y_hat = tf.nn.sigmoid(tf.add(tf.matmul(self.X,self.W),self.b))

        # Sigmoid Cross entropy cost function
        self.cost = tf.nn.sigmoid_cross_entropy_with_logits(logits = self.Y_hat, labels = self.Y)
        variable_summaries(self.cost,'Cross_Entropy_Cost')

        # Gradient Descent Optimizer
        self.optimizer = tf.train.GradientDescentOptimizer(learning_rate = self.alpha).minimize(self.cost)

        # Accuracy
        correct_prediction = tf.equal(tf.argmax(self.Y_hat,1),tf.argmax(self.Y,1))
        self.accuracy = tf.reduce_mean(tf.cast(correct_prediction,tf.float32))
        tf.summary.scalar('Accuracy',self.accuracy)

        # Global Variable Initiliazer
        self.init = tf.global_variables_initializer()

    def train(self):
        # starting tensorflow session
        with tf.Session() as sess:
            # Merge all the summaries and write them out to /tmp/mnist_logs (by default)
            merged = tf.summary.merge_all()
            train_writer = tf.summary.FileWriter('./train', sess.graph)
            test_writer = tf.summary.FileWriter('./test')


            # Initialize the variables
            sess.run(self.init)

            # create lists for storing cost and accuracy for every epoch
            cost_history, accuracy_history = [], []

            # iterate through all epochs
            for epoch in range(self.epochs):
                cost_per_epoch = 0

                # run the optimizer
                sess.run(self.optimizer, feed_dict = {self.X: self.X_train, self.Y: self.y_train})

                # calculate the cost of the current epoch
                c = sess.run(self.cost, feed_dict = {self.X: self.X_train, self.Y: self.y_train})

                # store the accuarcy and cost to the cost_history
                cost_history.append(sum(sum(c)))
                accuracy_history.append(self.accuracy.eval({self.X:self.X_train,self.Y:self.y_train})*100)

                # Displaying result on current Epoch
                if epoch % 25 == 0 and epoch != 0:
                    print("Epoch " + str(epoch) + " Cost: "
                                    + str(cost_history[-1]))

                # add sumaries to filewrite for tensorboard
                summary = sess.run(merged, feed_dict = {self.X:self.X_train,self.Y:self.y_train})
                train_writer.add_summary(summary,epoch)


            Weight = sess.run(self.W) # Optimized Weight
            Bias = sess.run(self.b)   # Optimized Bias

            # Final Accuracy
            print("\n Training Accuracy:", self.accuracy.eval({self.X:self.X_train,self.Y:self.y_train})*100, "%")
            #print("\n Validation Accuracy:", accuracy.eval({self.X:self.X_val,self.Y:self.y_val})*100, "%")
            print("\n Test Accuracy:", self.accuracy.eval({self.X:self.X_test,self.Y:self.y_test})*100, "%")
