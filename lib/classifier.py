"""implementation of a basic classifier class which allows to use different algorithms in the same framework."""
# classifiers
from sklearn.linear_model import LogisticRegression

# visulaization
from sklearn.metrics import classification_report

classifier_dict = {'lr':LogisticRegression}

class Classifier():
    def __init__(self, classifier='lr',**kwargs):
        """Initialization of a new Classifier object. This requires the name of classifier."""
        self.clf = classifier_dict[classifier](random_state=42,**kwargs)
        self.classes = ['non-pathogenic','pathogenic']
        self.trained = False

    def train(self,train_set):
        self.clf = self.clf.fit(train_set[0],train_set[1])
        print(train_set[1])
        y_pred = self.clf.predict(train_set[0])
        print(classification_report(train_set[1],y_pred,labels=self.classes))
        self.trained = True

    def test(self,test_set):
        if self.trained:
            y_pred = self.clf.predict(test_set[0])
            print(classification_report(test_set[1],y_pred,labels=self.classes))

    def fit(self,train_set,test_set):
        self.train(train_set)
        self.test(test_set)
