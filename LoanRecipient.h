//
// Created by Sam Clear on 10/21/21.
//

#ifndef CS124PROJECT4_LOANRECIPIENT_H
#define CS124PROJECT4_LOANRECIPIENT_H

#include <assert.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using std::string;
using std::ostream;
using std::iostream;
using std::ifstream;
using std::vector;
using std::setw;
using std::right;
using std::left;

using namespace std;
/**
 * This class presents data compiled on Predictions of clients defaulting on loans based on client behaviors, with the
 * data generated from Univ.AI for the Geoffrey Hinton Scholarship's Hackathon competition. Information on the event and
 * the organization running this competition can be found at https://www.univ.ai/ghf, and information for the data can
 * be found at the Univ.AI data set description here https://hack.univ.ai/problem/data.
 *
 * Fields (names abbreviated as represented in the fields below)
 *  id: An integer serial number representing the loan recipient's place in the list
 *  income: An integer representing the total income of the recipient
 *  age: An integer representing the age of the recipient
 *  YoE (Years of Job Experience): An integer representing the professional experience of the user in years
 *  mar: A string representing the marital status of the recipient
 *  houseOwn: A string representing the house ownership status of the recipient
 *  carOwn: A string representing the car ownership status of the recipient
 *  job: A string representing the current profession of the recipient
 *  city: A string representing the current city of the recipient's residence. In this field variable, specific cities
 *      that are two words have the space replaced with an underscore, and cities with commas in their name have the
 *      commas removed.
 *  state: A string representing the current state of the recipient's residence
 *  currJobYrs: An integer variable representing the number of years the recipient has held their current job
 *  currHouseYrs: An integer variable representing the number of years the recipient has owned their current house
 *  riskFlag: A boolean variable that represents whether the recipient has defaulted on a loan or not
 */

class LoanRecipient {
    // Fields
private:
    int id, income, age, YoE, currJobYrs, currHouseYrs;
    string mar, houseOwn, carOwn, job, city, state;
    int riskFlag;

public:
    // Default Constructor Override
    LoanRecipient() {
        id = 0;
        income = 0;
        age = 0;
        YoE = 0;
        mar = "";
        houseOwn = "";
        carOwn = "";
        job = "";
        city = "";
        state = "";
        currJobYrs = 0;
        currHouseYrs = 0;
        riskFlag = 0;
    }

    // Constructor Method
    LoanRecipient(int id, int income, int age, int YoE, string mar, string houseOwn, string gCarOwn, string job,
                  string city, string state, int currJobYrs, int currHouseYrs, bool riskFlag) {
        this->id = id;
        this->income = income;
        this->age = age;
        this->YoE = YoE;
        this->mar = mar;
        this->houseOwn = houseOwn;
        carOwn = gCarOwn;
        this->job = job;
        this->city = city;
        this->state = state;
        this->currJobYrs = currJobYrs;
        this->currHouseYrs = currHouseYrs;
        this->riskFlag = riskFlag;
    }

    // Getter / Accessor Methods
    int getID() const {
        return id;
    }

    int getIncome() const {
        return income;
    }

    int getAge() const {
        return age;
    }

    int getYearsOfExp() const {
        return YoE;
    }

    string getMarStatus() const {
        return mar;
    }

    string getHouseOwn() const {
        return houseOwn;
    }

    string getCarOwn() const {
        return carOwn;
    }

    string getJob() const {
        return job;
    }

    string getCity() const {
        return city;
    }

    string getState() const {
        return state;
    }

    int getCurrJobYrs() const {
        return currJobYrs;
    }

    int getCurrHouseYrs() const {
        return currHouseYrs;
    }

    bool getRiskFlag() const {
        return riskFlag;
    }

    // Setter / Mutator Methods
    void setID(int id) {
        this->id = id;
    }

    void setIncome(int income) {
        this->income = income;
    }

    void setAge(int age) {
        this->age = age;
    }

    void setYearsOfExp(int YoE) {
        this->YoE = YoE;
    }

    void setMarStatus(string mar) {
        this->mar = mar;
    }

    void setHouseOwn(string houseOwn) {
        this->houseOwn = houseOwn;
    }

    void setCarOwn(string carOwn) {
        this->carOwn = carOwn;
    }

    void setJob(string job) {
        this->job = job;
    }

    void setCity(string city) {
        this->city = city;
    }

    void setState(string state) {
        this->state = state;
    }

    void setCurrJobYears(int currJobYrs) {
        this->currJobYrs = currJobYrs;
    }

    void setCurrHouseYrs(int currHouseYrs) {
        this->currHouseYrs = currHouseYrs;
    }

    void setRiskFlag(bool riskFlag) {
        this->riskFlag = riskFlag;
    }

    // This code was taken from the Project 1 Example given to us under the Project 1 Folder in BlackBoard. This segment
    // is authored by Professor Lisa Dion and is altered to fit the LoanRecipient class.
    /**
     * This overloaded "<" operator allows us to compare LoanRecipient objects when they must be evaluated for
     * object comparison and placement on the Binary Search Tree
     * @param lhs LoanRecipient object on the left hand side of the evaluation
     * @param rhs LonaRecipient object on the right hand side of the evaluation
     */
    friend bool operator < (const LoanRecipient &lhs, const LoanRecipient &rhs) {
        return lhs.id < rhs.id;
    }

    /**
     * This overloaded "<" operator allows us to compare LoanRecipient objects when they must be evaluated for
     * object comparison and placement on the Binary Search Tree
     * @param lhs LoanRecipient object on the left hand side of the evaluation
     * @param rhs LonaRecipient object on the right hand side of the evaluation
     */
    friend bool operator > (const LoanRecipient &lhs, const LoanRecipient &rhs) {
        return lhs.id > rhs.id;
    }

    /**
     * This overloaded "<" operator allows us to compare LoanRecipient objects when they must be evaluated for
     * object comparison and placement on the Binary Search Tree
     * @param lhs LoanRecipient object on the left hand side of the evaluation
     * @param rhs LonaRecipient object on the right hand side of the evaluation
     */
    friend bool operator <= (const LoanRecipient &lhs, const LoanRecipient &rhs) {
        return lhs.id <= rhs.id;
    }

    /**
     * This overloaded "<" operator allows us to compare LoanRecipient objects when they must be evaluated for
     * object comparison and placement on the Binary Search Tree
     * @param lhs LoanRecipient object on the left hand side of the evaluation
     * @param rhs LonaRecipient object on the right hand side of the evaluation
     */
    friend bool operator >= (const LoanRecipient &lhs, const LoanRecipient &rhs) {
        return lhs.id >= rhs.id;
    }

    /**
     * This overloaded "==" operator allows us to compare LoanRecipient objects when they must be evaluated for
     * object comparison and placement on the Binary Search Tree
     * @param lhs LoanRecipient object on the left hand side of the evaluation
     * @param rhs LonaRecipient object on the right hand side of the evaluation
     */
    friend bool operator == (const LoanRecipient &lhs, const LoanRecipient &rhs) {
        // Comparing two objects unique Income variables
        return lhs.id == rhs.id;
    }

    /**
     * This function checks the data set's fields to determine whether there are 252,000 Loan Recipient objects with
     * valid field values.
     * @param loaners The vector of LoanRecipient object being tested
     * @return Boolean variable indicating whether the vector contains 252,000 sets of valid fields
     */
    bool fieldValidation(vector<LoanRecipient> loaners) {
        // First, vector is checked to have the correct number of rows for each of the loan recipients present in the
        // data set. If there are less than the total number of loan recipients, minus 1 to account for index starting
        // at zero, then the function immediately returns false as to avoid crashing of the program.
        if (loaners.size() < 251999) {
            return false;
        }

        // Secondly, index starts at zero and the dataStatus variable is set to true. Then, a while loop runs through
        // each loan recipient object to check for an invalid field value until it reaches the end of the vector. If an
        // invalid field value is found, then the data is assumed invalid and the method returns false.
        int idx = 0;
        bool dataStatus = true;
        while(idx < loaners.size() && dataStatus == true) {

            // Check that each index value is one less than the object's ID, return false and end loop if not.
            if (loaners[idx].getID() != idx + 1) {
                dataStatus = false;
            }

            // Check that each income value is positive, return false and end loop if not.
            if (loaners[idx].getIncome() < 0) {
                dataStatus = false;
            }

            // Check that each age value is positive, return false and end loop if not.
            if (loaners[idx].getAge() < 0) {
                dataStatus = false;
            }

            // Check that each Years of Experience value is positive, return false and end loop if not.
            if (loaners[idx].getYearsOfExp() < 0) {
                dataStatus = false;
            }

            // Check that each Marital Status value is either "single" or "married", return false and end loop if not.
            if (loaners[idx].getMarStatus().compare("single") != 0 &&
                loaners[idx].getMarStatus().compare("married") != 0) {
                dataStatus = false;
            }

            // Check that each House Ownership entry is either "rented", "owned", or "norent_noown", return false and
            // end loop if not.
            if (loaners[idx].getHouseOwn().compare("rented") != 0 &&
                loaners[idx].getHouseOwn().compare("owned") != 0 &&
                loaners[idx].getHouseOwn().compare("norent_noown") != 0) {
                dataStatus = false;
            }

            // Check that each Car Ownership entry is either "yes" or "no", return false and end loop if not.
            if (loaners[idx].getCarOwn().compare("yes") != 0 && loaners[idx].getCarOwn().compare("no") != 0) {
                dataStatus = false;
            }

            // Check that the Profession of the LoanRecipient object is not empty, return false and end loop if not.
            if (loaners[idx].getJob().compare("") == 0) {
                dataStatus = false;
            }

            // Check that the City of the LoanRecipient object is not empty, return false and end loop if not.
            if (loaners[idx].getCity().compare("") == 0) {
                dataStatus = false;
            }

            // Check that the State of the LoanRecipient object is not empty, return false and end loop if not.
            if (loaners[idx].getState().compare("") == 0) {
                dataStatus = false;
            }

            // Check that each Current Job Years value is positive, return false and end loop if not.
            if (loaners[idx].getCurrJobYrs() < 0) {
                dataStatus = false;
            }

            // Check that every Current House Years value is positive, return false and end loop if not.
            if (loaners[idx].getCurrHouseYrs() < 0) {
                dataStatus = false;
            }

            // Check that each Risk Flag boolean value is either true or false, return false and end loop if not.
            if (loaners[idx].getRiskFlag() != true && loaners[idx].getRiskFlag() != false) {
                dataStatus = false;
            }
            ++idx;
        }

        // Finally, the status of the car ownership data is return to the method call.
        return dataStatus;
    }

    /**
     * This overloaded ostream operator acts as the printing method for each LoanRecipient object. This overload of the
     * ostream returns a string of each of the values present in the LoanRecipient object with proper spacing to make a
     * table.
     * @param outs A ostream object containing the different fields, spacings, and attributes of the LoanRecipient
     *             object
     * @param loaned A reference to the LoanRecipient object being printed
     * @return An ostream object representing the fields in the LoanRecipient object in the given format
     */
    friend ostream& operator << (ostream& outs, const LoanRecipient &loaned) {
        outs << left << setw(5)<< loaned.id;
        /**
        outs << left << setw(10) << loaned.income;
        outs << setw(5) << loaned.age;
        outs << setw(5) << loaned.YoE;
        outs << setw(15) << loaned.mar;
        outs << setw(10) << loaned.houseOwn;
        outs << setw(5) << loaned.carOwn;
        outs << setw(20) << loaned.job;
        outs << setw(20) << loaned.city;
        outs << setw(20) << loaned.state;
        outs << setw(6) << loaned.currJobYrs;
        outs << setw(8) << loaned.currHouseYrs;
        outs << setw(5) << loaned.riskFlag;
         */
        return outs;
    }

/**
  * This method tests the accessor and mutator methods of the LoanRecipient object class.
  * @return Boolean value indicating success of methods tested
  */
    bool loanMethodTester() {
        // Create a LoanRecipient object to test accessor and mutator methods on using constructor
        LoanRecipient tester = LoanRecipient(1,1,1,1,"single","rented",
                                             "no","Software_engineer","Rewa","Madhya_Pradesh",
                                             1,1,1);
        // Each accessor method is tested with the given tester fields
        assert(tester.getID() == 1);
        assert(tester.getIncome() == 1);
        assert(tester.getAge() == 1);
        assert(tester.getYearsOfExp() == 1);
        assert(tester.getMarStatus().compare("single") == 0);
        assert(tester.getHouseOwn().compare("rented") == 0);
        assert(tester.getCarOwn().compare("no") == 0);
        assert(tester.getCity().compare("Rewa") == 0);
        assert(tester.getState().compare("Madhya_Pradesh") == 0);
        assert(tester.getCurrJobYrs() == 1);
        assert(tester.getCurrHouseYrs() == 1);
        assert(tester.getRiskFlag() == 1);

        // Then, each field is changed and the accessor methods evaluate whether the setter methods were successful
        tester.setID(2);
        assert(tester.getID() == 2);
        tester.setIncome(2);
        assert(tester.getIncome() == 2);
        tester.setAge(2);
        assert(tester.getAge() == 2);
        tester.setYearsOfExp(2);
        assert(tester.getYearsOfExp() == 2);
        tester.setMarStatus("married");
        assert(tester.getMarStatus().compare("married") == 0);
        tester.setHouseOwn("owned");
        assert(tester.getHouseOwn().compare("owned") == 0);
        tester.setCarOwn("yes");
        assert(tester.getCarOwn().compare("yes") == 0);
        tester.setCity("New York City");
        assert(tester.getCity().compare("New York City") == 0);
        tester.setState("New York");
        assert(tester.getState().compare("New York") == 0);
        tester.setCurrJobYears(2);
        assert(tester.getCurrJobYrs() == 2);
        tester.setCurrHouseYrs(2);
        assert(tester.getCurrHouseYrs() == 2);
        tester.setRiskFlag(false);
        assert(tester.getRiskFlag() == 0);
        return true;
    }
};

/** This function calculates the mean of all the Loan Recipient's bank accounts. This is done by adding each of the
     *  Loaner Recipient object's income values into the accumulator field incomeTotal, and dividing the total income
     *  accumulated by the total number of LoanRecipient objects. Initially, the int datatype could not contain the
     *  sum of each of the income totals, which lead me to use the long datatype to contain the entirety of the total
     *  income of all the LoanRecipient objects. The long datatype was pulled from my previous experience with the
     *  different int data types in CS 110.
     * @param loanerList Vector containing each of the LoanerRecipient object's locations
     * @return Average income of each Loan Recipient in the vector given
     */
double getIncomeMean(vector<LoanRecipient> loaners) {
    long incomeTotal = 0, indexMax = loaners.size(), count = 0;
    while (count < indexMax) {
        incomeTotal = incomeTotal + loaners[count].getIncome();
        ++count;
    }
    return incomeTotal / indexMax;
}

/**
 * This method initializes a given LoanRecipient vector with initialized LoanRecipient objects.
 * @param filename A string indicating the location of the file of LoanRecipient field values
 * @param loaners An empty vector with a LoanRecipient datatype given from the main function
 * @return A boolean value indicating whether the file was successfully accessed and the vector was successfully
 *         initialized with LoanRecipient objects.
 */

bool getDataFromFile(string filename, vector<LoanRecipient> &loaners) {
    bool ran;

    ifstream FIn;
    FIn.open("../Training Data.csv");

    // Load in headers into Headers Var
    string headers;
    if (FIn) {
        getline(FIn, headers);
        // std::cout << headers << std::endl;
    }

    // Begin loading data into the vector
    while (FIn && FIn.peek() != EOF) {

        // Each field for each Loan Recipient
        int id, income, age, YoE, currJobYrs, currHouseYrs;
        string mar, houseOwn, carOwn, job, city, state;
        string city1, city2; // These fields are used to connect city names seperated by commas
        bool riskFlag;
        char comma, quote; // This field is used as a placeholder to divide variable values by commas and quotes in the list

        // Read in ID field, Differentiate Variables by Commas
        FIn >> id;
        FIn >> comma;

        // Read in income field, then comma
        FIn >> income;
        FIn >> comma;

        // Read in age field, then comma
        FIn >> age;
        FIn >> comma;

        // Read in Years of Experience field, then comma
        FIn >> YoE;
        FIn >> comma;

        // Read in Marriage Status field, then comma
        getline(FIn, mar, ',');

        // Read in House Ownership Status field, then comma
        getline(FIn, houseOwn, ',');

        // Read in Car Ownership Status field, then comma
        getline(FIn, carOwn, ',');

        // Read in Job Title field, then comma
        getline(FIn, job, ',');

        // Read in City field, then comma
        if (FIn.peek() == '"') {
            FIn >> quote;
            getline(FIn, city1, ',');
            FIn >> comma;
            getline(FIn, city2, ',');
            city = city1 + city2;
        } else {
            getline(FIn, city, ',');
        }

        // Read in State field, then comma
        getline(FIn, state, ',');

        // Read in Years in Current Job field, then comma
        FIn >> currJobYrs;
        FIn >> comma;

        // Read in Years of House Ownership field, then comma
        FIn >> currHouseYrs;
        FIn >> comma;

        // Read in Risk Flag status, then comma
        FIn >> riskFlag;

        // Finally, read in the fields into a Loan Recipient object in the vector
        loaners.push_back(LoanRecipient(id, income, age, YoE, mar, houseOwn, carOwn, job, city, state, currJobYrs, currHouseYrs,riskFlag));
    }
    ran = true;
    return ran;
}

#endif //CS124PROJECT4_LOANRECIPIENT_H
