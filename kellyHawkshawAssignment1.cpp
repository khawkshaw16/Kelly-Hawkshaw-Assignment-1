#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>
#include <math.h>
using namespace std;
//set up the bigram calc for the nucleotides

int main(int argc, char **argv){
  //used for user if they want to execute the program continuously
  bool executeCont = true;
  //asking for file
  while(executeCont) {
    //asking user for the text file name as input
    string nameOfFile;
    cout << "Enter name of File: " << endl;
    cin >> nameOfFile;
    ifstream inTextFile;
    inTextFile.open(nameOfFile);

    //Checking to make sure the file opens properly before proceeding
    if(inTextFile.fail()){

      cerr << "Text file opening error!!!" << endl;
      exit(1);
    }

    //Setting up DNA string variables
    string dnaString;
    float sumDNA;
    float count;
    float mean;
    float variance;
    float stdDev;
    float max;
    float maxLength;
    float tempMax;
    float tempSum;

    int numA;
    int numC;
    int numG;
    int numT;

    float probabilityOfA;
    float probabilityOfC;
    float probabilityOfG;
    float probabilityOfT;

    float probabilityOfAA;
    float probabilityOfAC;
    float probabilityOfAG;
    float probabilityOfAT;
    float probabilityOfCA;
    float probabilityOfCC;
    float probabilityOfCG;
    float probabilityOfCT;
    float probabilityOfGA;
    float probabilityOfGC;
    float probabilityOfGG;
    float probabilityOfGT;
    float probabilityOfTA;
    float probabilityOfTC;
    float probabilityOfTG;
    float probabilityOfTT;

    //Reading in the values from text file
      while(!inTextFile.eof()){

        inTextFile >> dnaString;
        tempMax = dnaString.length();
        //check maxLength
        if(tempMax > maxLength){

          maxLength = tempMax;
        }

        sumDNA += dnaString.length();
        count++;

        //getting mean
        mean = sumDNA/count;

        //getting the tempSum
        tempSum += pow((dnaString.length() - mean), 2);

        //Getting the number of each of the letters (ACGT)
        for(int i = 0; i < dnaString.length(); ++i){

          if(dnaString[i] == 'A'){

            numA++;
          }

          if(dnaString[i] == 'C'){

            numC++;
          }

          if(dnaString[i] == 'G'){

            numG++;
          }

          if(dnaString[i] == 'T'){

            numT++;
          }

        }

      }

      //Getting the probablity of the letters (ACGT) and then the nucleotides
      probabilityOfA = numA/sumDNA;
      probabilityOfC = numC/sumDNA;
      probabilityOfG = numG/sumDNA;
      probabilityOfT = numT/sumDNA;
      probabilityOfAA = probabilityOfA*probabilityOfA;
      probabilityOfAC = probabilityOfA*probabilityOfC;
      probabilityOfAG = probabilityOfA*probabilityOfG;
      probabilityOfAT = probabilityOfA*probabilityOfT;
      probabilityOfCA = probabilityOfC*probabilityOfA;
      probabilityOfCC = probabilityOfC*probabilityOfC;
      probabilityOfCG = probabilityOfC*probabilityOfG;
      probabilityOfCT = probabilityOfC*probabilityOfT;
      probabilityOfGA = probabilityOfG*probabilityOfA;
      probabilityOfGC = probabilityOfG*probabilityOfC;
      probabilityOfGG = probabilityOfG*probabilityOfG;
      probabilityOfGT = probabilityOfG*probabilityOfT;
      probabilityOfTA = probabilityOfT*probabilityOfA;
      probabilityOfTC = probabilityOfT*probabilityOfC;
      probabilityOfTG = probabilityOfT*probabilityOfG;
      probabilityOfTT = probabilityOfT*probabilityOfT;

      //Getting standard deviation
      stdDev = sqrt(variance);

      //Getting the variance
      variance = tempSum/count - 1;
      variance *= -1;

      //outputing the text file
      ofstream outputTextFile;
      outputTextFile.open("kellyhawkshaw.txt");
      outputTextFile << "Kelly Hawkshaw ID: 2328274" << endl;
      outputTextFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
      outputTextFile << "The sum legth of the DNA string: " << sumDNA << endl;
      outputTextFile << "The mean: " << mean << endl;
      outputTextFile << "The standard deviation: " << stdDev <<endl;
      outputTextFile << "The variance: " << variance << endl;
      outputTextFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
      outputTextFile << "Probability of A: " << probabilityOfA << endl;
      outputTextFile << "Probability of C: " << probabilityOfC << endl;
      outputTextFile << "Probability of G: " << probabilityOfG << endl;
      outputTextFile << "Probability of T: " << probabilityOfT << endl;
      outputTextFile << "Probability of AA: " << probabilityOfAA << endl;
      outputTextFile << "Probability of AC: " << probabilityOfAC << endl;
      outputTextFile << "Probability of AG: " << probabilityOfAG << endl;
      outputTextFile << "Probability of AT: " << probabilityOfAT << endl;
      outputTextFile << "Probability of CA: " << probabilityOfCA << endl;
      outputTextFile << "Probability of CC: " << probabilityOfCC << endl;
      outputTextFile << "Probability of CG: " << probabilityOfCG << endl;
      outputTextFile << "Probability of CT: " << probabilityOfCT << endl;
      outputTextFile << "Probability of GA: " << probabilityOfGA << endl;
      outputTextFile << "Probability of GC: " << probabilityOfGC << endl;
      outputTextFile << "Probability of GG: " << probabilityOfGG << endl;
      outputTextFile << "Probability of GT: " << probabilityOfGT << endl;
      outputTextFile << "Probability of TA: " << probabilityOfTA << endl;
      outputTextFile << "Probability of TC: " << probabilityOfTC << endl;
      outputTextFile << "Probability of TG: " << probabilityOfTG << endl;
      outputTextFile << "Probability of TT: " << probabilityOfTT << endl;
      outputTextFile << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
      outputTextFile << "1000 DNA strings using Gaussian distribution: " <<endl;

      string gdDNAString;

      int nucleotideA;
      int nucleotideC;
      int nucleotideG;
      int nucleotideT;

      //Iterating the 100 dna string using gaussian distribution
      for(int i = 0; i < 1000; ++i){
        //gaussian distribution
        float randA = (float) rand()/RAND_MAX;
        float randB = (float) rand()/RAND_MAX;
        float gdEquation = (sqrt(-2 * log(randA)) * cos(2 * M_PI * randB));
        int gd = stdDev * gdEquation + mean;

        float nucleotideA = numA/sumDNA * gd;
        float nucleotideC = numC/sumDNA * gd;
        float nucleotideG = numG/sumDNA * gd;
        float nucleotideT = numT/sumDNA * gd;

        //Making the string iteration for the letters
        for(int n = 0; n < nucleotideA; ++n){

          gdDNAString += "A";
        }
        for(int n = 0; n < nucleotideC; ++n){

          gdDNAString += "C";
        }
        for(int n = 0; n < nucleotideG; ++n){

          gdDNAString += "G";
        }
        for(int n = 0; n < nucleotideT; ++n){

          gdDNAString += "T";
        }

        outputTextFile << gdDNAString << endl;

      }

      //Closing the files then asking if they would like to do it again
      inTextFile.close();
      outputTextFile.close();
      cout << "Would you like to process another list, if so type Y for yes and N for no" << endl;
      //if loop for executing the program again
      char userInput;
      if(userInput == 'N' || 'n') {

        executeCont = true;
      }
      else if(userInput == 'N' || 'n'){

        executeCont = false;
        cout << "You have exited the program have a good day!" << endl;
        break;
      }
      else {

        executeCont = false;
        cout << "Error occured exiting program have a good day!" << endl;
        break;
      }

  }

  return 0;

}
