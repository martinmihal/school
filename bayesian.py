#author: Martin Mihal
#about: python program for detect language
import os
from collections import  defaultdict
import random
import numpy as np
import  sys

models = [] # list for storing all models

class AllLanguageData():
    """
    class which store information about all languages and is able to predict language
    """
    languagesDataObjects = dict()

    def __init__(self):
        self.languagesDataObjects =  dict()

    def addLanguage(self, languageDataObject):
        self.languagesDataObjects[languageDataObject.nameOfLanguage] = languageDataObject

    def getLanguage(self, nameOfLanguage):
        return self.languagesDataObjects.get(nameOfLanguage)

    def getNamesOfLanuages(self):
        return self.languagesDataObjects.keys()

    def decideLanguage(self, listOfWords):
        """
        method which count language probalities for given sentence
        and winning(max) probality is returned as recognized language
        """
        maxProbalityYet, winnerLanguage = -1, None
        for language in self.getNamesOfLanuages():
            probForLanguage = 1
            languageData = self.languagesDataObjects.get(language)
            for word in listOfWords:
                probForLanguage = probForLanguage * languageData.laplacianSmothingProbability(word)
            percentageOfLanguageInTraining =  self.getNumberOfAllFilesForLanguage(language) / float(self.getNumberOfAllFiles())
            probForLanguage = probForLanguage * percentageOfLanguageInTraining
            if probForLanguage > maxProbalityYet:
                maxProbalityYet = probForLanguage
                winnerLanguage = language
        return winnerLanguage

    def getNumberOfAllFilesForLanguage(self, language):
        return  self.languagesDataObjects.get(language).numberOfFiles

    def getNumberOfAllFiles(self):
        result = 0
        for language in self.getNamesOfLanuages():
            languageData = self.languagesDataObjects.get(language)
            result += languagesData.numberOfFiles
        return result

    def testAccurancyOnTrainData(self):
        '''
        :return:ratio of true and false recognized items on train data
        '''
        trueRecognized, falseRecognized = 0, 0
        for language in self.getNamesOfLanuages():
            languageData = self.getLanguage(language)
            for file in languageData.trainData:
                pathToFile =  os.path.join(languageData.folderForLanguage, file)
                with open(pathToFile) as f:
                    words = open(pathToFile).read().split()
                    recognizedLanguage = database.decideLanguage(words)
                    if language == recognizedLanguage:
                        trueRecognized += 1
                    else:
                        falseRecognized += 1
        return trueRecognized / float(trueRecognized + falseRecognized)

    def testAccurancyOnTestData(self):
        '''
        :return:ratio of true and false recognized items on test data
        '''
        trueRecognized, falseRecognized = 0, 0
        for language in self.getNamesOfLanuages():
            languageData = self.getLanguage(language)
            for file in languageData.testData:
                pathToFile =  os.path.join(languageData.folderForLanguage, file)
                with open(pathToFile) as f:
                    words = open(pathToFile).read().split()
                    recognizedLanguage = database.decideLanguage(words)
                    if language == recognizedLanguage:
                        trueRecognized += 1
                    else:
                        falseRecognized += 1
        return trueRecognized / float(trueRecognized + falseRecognized)


class LanguageData():
    """
    class which store learned data for one language
    """
    wordInLanguage = defaultdict(int)  #count number of occurence for every word in given language

    numberOfFiles = 0

    nameOfLanguage = None

    trainData = []

    testData = []

    folderForLanguage = None

    def __init__(self, nameOfLanguage, numberOfFiles, wordInLanguage, trainData, testData, folderForLanguage):
        self.nameOfLanguage = nameOfLanguage
        self.numberOfFiles = numberOfFiles
        self.wordInLanguage = wordInLanguage
        self.trainData = trainData
        self.testData = testData
        self.folderForLanguage = folderForLanguage

    def getNumberOfAllUsedWords(self):
        return sum(self.wordInLanguage.values())

    def laplacianSmothingProbability(self, word):
        """
            :return probability that give word is from give language
        """
        numberOfOccurenceInCorpus = self.wordInLanguage.get(word)
        if numberOfOccurenceInCorpus is None:
            numberOfOccurenceInCorpus = 0
        allSeenFeatures = len(self.wordInLanguage)
        return numberOfOccurenceInCorpus + 1 / float(allSeenFeatures + self.numberOfFiles)

def getTrainAndTestData(folderForLanguage, whichPart):
    allFiles = os.listdir(folderForLanguage)
    howManyFilesInOnePart = 10
    splittedFiles = np.array_split(np.array(allFiles), howManyFilesInOnePart) # split data into howManyFilesInOnePart parts
    trainData, testData = [], []
    for i, subList in enumerate(splittedFiles):
        if i == whichPart:          # one part of splited list is used as test data
            testData.extend(subList)
        else:                       # all others as train data
            trainData.extend(subList)
    return trainData, testData

def cleanInput(listOfWords):
    """
    for better learning and recognition we want ignore these special characaters
    :param listOfWords:
    :return cleared words
    """
    result = []
    for word in listOfWords:
        correctWord = word.replace(":",'')
        correctWord = correctWord.replace(".",'')
        correctWord = correctWord.replace(",",'')
        correctWord = correctWord.replace("!",'')
        result.append(correctWord)
    return  result

#loading of data from file
currentPath =  os.getcwd()
dataFolder = os.path.join(currentPath, "data")
languages = os.listdir(dataFolder)

########### learning of models ###########
for kFoldVAlue in range(10):
    database = AllLanguageData()
    for language in languages:
        wordInLanguage = defaultdict(int)  #count for every word in given language
        numberOfFileForLanguage = 0
        folderForLanguage = os.path.join(dataFolder, language)
        trainData, testData = getTrainAndTestData(folderForLanguage, kFoldVAlue) # divide data
        for file in trainData:
            numberOfFileForLanguage += 1
            pathToFile =  os.path.join(folderForLanguage, file)
            with open(pathToFile) as f:
                words = open(pathToFile).read().split()
                words = cleanInput(words)
                for word in words:
                    wordInLanguage[word]+= 1
        languagesData = LanguageData(language, numberOfFileForLanguage, wordInLanguage,  trainData, testData, folderForLanguage)
        database.addLanguage(languagesData) #store all learned data about language into main object
    models.append(database) #added this model with given k-fold value into set of models

########### selecting best model - that one with max accurancy on test data ###########
maxAccurancyYet, winningModel = 0, None
for model in models:
    accurancyOnTestData = model.testAccurancyOnTestData()
    print(accurancyOnTestData)
    if accurancyOnTestData > maxAccurancyYet:
        maxAccurancyYet = accurancyOnTestData
        winningModel = model

print("accurancy of winning model on train data", winningModel.testAccurancyOnTrainData()) # it's part of score for this project
########### printing of results ###########
while 1:
    try:
        line = sys.stdin.readline()
    except KeyboardInterrupt:
        break
    if not line:
        break
    rowToListOfWords = cleanInput(line.split(" "))
    recognizedLanguage = winningModel.decideLanguage(rowToListOfWords)
    print(recognizedLanguage)