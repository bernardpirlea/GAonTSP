#include <iostream>
#include <time.h>
#include <math.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <random>

using namespace std;

int randomInRange(int range_from,int range_to)
{
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(range_from, range_to);
    return distr(generator);
}

// Function to calculate distance
int distance2D(int x1, int y1, int x2, int y2)
{
    // Calculating distance
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

void displaySolution(vector<int> vect){
    for (int i = 0; i < vect.size(); i++) {
        cout<< vect.at(i)<< " ";
    }
    cout<<endl;
}
void displayFitness(vector<float> vect){
    for (int i = 0; i < vect.size(); i++) {
        cout<< vect.at(i)<<endl;
    }
}
void displayPopulation(vector<vector<int> > population)
{
    for(int i = 0; i < population.size(); i++)
    {
        displaySolution(population.at(i));
    }
    cout<<endl;
}

int evaluate_solution(vector<vector<float >> coord, vector<int> solution)
{
    vector<float> point1;
    vector<float> point2;
    int distance = 0;
    for(int j = 0; j < solution.size() - 1; j++)
    {
        point1 = coord.at(solution.at(j));
        point2 = coord.at(solution.at(j + 1));
        distance += distance2D(point1.at(0), point1.at(1), point2.at(0), point2.at(1));
        //cout<< point1.at(0) << " " << point1.at(1)<<endl;
        //cout<< point2.at(0) << " " << point2.at(1)<<endl;
        //cout<< distance<<endl;
    }
    point1 = coord.at(solution.at(solution.size() - 1));
    point2 = coord.at(solution.at(0));
    distance += distance2D(point1.at(0), point1.at(1), point2.at(0), point2.at(1));
    
    return distance;
}

void readData(string fileName, vector<vector<float >> & coord)
{
    string delimiter = " ";
    string line;
    ifstream myfile (fileName);
    if (myfile.is_open())
    {
        while (getline (myfile, line))
        {
            vector<float> point;
            int i = 0;
            size_t pos = 0;
            string token;
            while ((pos = line.find(delimiter)) != string::npos)
            {
                token = line.substr(0, pos);
                if(i > 0)
                {
                    point.push_back(stof(token));
                }
                line.erase(0, pos + delimiter.length());
                i++;
            }
            point.push_back(stof(line));
            //cout << line << '\n';
            coord.push_back(point);
        }
        myfile.close();
    }
    else cout << "Unable to open file";
}

vector<int> generate_solution(unsigned long nrOfCities)
{
    vector<int> solution;
    std::random_device rd;
    std::default_random_engine engine{rd()};
    for(int i = 0; i < nrOfCities; i++)
    {
        solution.push_back(i);
    }
    
    shuffle(std::begin(solution), std::end(solution), engine);
    return solution;
}

vector<vector<int> > generate_population(int pop_size, unsigned long nrOfCities)
{
    std::random_device rd;
    std::default_random_engine engine{rd()};
    vector<vector<int> > population;
    vector<int> order;
    
    for(int i = 0; i < nrOfCities; i++)
    {
        order.push_back(i);
    }
    
    population.push_back(order);
    
    for(int p = 1; p < pop_size; ++p)
    {
        //randomly generate one individual
        std::shuffle(std::begin(order), std::end(order), engine);
        
        //add the individual to the population
        population.push_back(order);
    }
    return population;
}


float random01()
{
    return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}

int randomRange(unsigned long max)
{
    return rand() % max;
}

vector<float> calculate_fitness(vector<vector<int> > population, vector<vector<float >> coord)
{
    vector<float> fitness;
    vector<float> point1;
    vector<float> point2;
    vector<float> functValue;
    float maxVal = - 1;
    
    //calculate max of population
    for( int i = 0; i < population.size(); i++)
    {
        int distance = 0;
        distance = evaluate_solution(coord, population.at(i));
        if(maxVal < distance)
        {
            maxVal = distance;
        }
        functValue.push_back(distance);
        //cout<< distance << endl;
    }
    
    //calculate fitness for each individual
    for(int i = 0; i < population.size(); ++i)
    {
        fitness.push_back(1.1 * maxVal - functValue.at(i));
    }

    return fitness;
}

vector<float> calculate_sf(vector<float> fitness)
{
    vector<float> sf;
    
    sf.push_back(fitness.at(0));
    for(int i = 1; i < fitness.size(); ++i)
    {
        //compute probability of selecting individual i
        sf.push_back(sf.at(i-1) + fitness.at(i));
    }
    return sf;
}

int selectF(vector<float> sf)
{
    float pos =  random01() * sf.at(sf.size()-1);
    
    int i;
    for(i = 0; i < sf.size(); ++i)
    {
        if(pos <= sf.at(i))
            return i;
    }
    return i;
}

vector<vector<int> > selection(vector<vector<int> > population, vector<float> fitness)
{
    vector<vector<int> > new_pop;
    vector<float> sf = calculate_sf(fitness);
    vector<int> index;
    for(int i = 0; i < population.size(); ++i)
    {
        int j = selectF(sf);
        new_pop.push_back(population.at(j));
    }
    return new_pop;
}

int get_position_best(vector<float> fitness)
{
    float maxValue = -1;
    int pos = 0;
    
    for(int i = 0; i < fitness.size(); ++i)
    {
        if(fitness.at(i) > maxValue)
        {
            maxValue = fitness.at(i);
            pos = i;
        }
    }
    return pos;
}

int get_best_value(vector<vector<int> > population, vector<float> fitness, vector<vector<float >> coord)
{
    vector<float> point1;
    vector<float> point2;
    vector<int> individual;
    int distance = 0;
    int posmin;
    posmin = get_position_best(fitness);
    individual = population.at(posmin);
    
    distance = evaluate_solution(coord, individual);
    
    return distance;
}

void mutation(vector<vector<int> >& population, int generation)
{
    vector<int> individual;
    float pm = 0.02;
    
    for(int i = 0; i < population.size(); ++i)
    {
        individual = population.at(i);
        for(int j = 0; j < individual.size(); ++j)
        {
            int point1 = randomRange(individual.size());
            int point2 = randomRange(individual.size());
            if(random01() < pm)
            {
                //swap the 2 points in the order
                int aux = individual.at(point1);
                individual.at(point1) = individual.at(point2);
                individual.at(point2) = aux;
                
            }
        }
        population.at(i) = individual;
        individual.clear();
    }
}

void crossover(vector<vector<int> >& population, int generation)
{
    float pc = 0.7;
    vector<int> pos;
    vector<int> individual1;
    vector<int> individual2;
    
    
    //choose random individuals
    for(int i = 0; i < population.size(); ++i)
    {
        if(random01() < pc)
        {
            pos.push_back(i);
        }
    }
    
    //make them pairs
    if(pos.size() % 2 == 1)
    {
        pos.pop_back();
    }
    
    
    //do cross over
    for(int i = 0; i < pos.size(); i +=2)
    {
        vector<int> newIndividual1;
        vector<int> newIndividual2;
        
        individual1 = population.at(pos.at(i));
        individual2 = population.at(pos.at(i+1));
        
        int slice = rand() % (individual1.size() - 2) + 1;
        
        
        //take first half of both individual and put it in the new individual1
        for(int j = 0; j < slice; j++)
        {
            newIndividual1.push_back(individual1.at(j));
            newIndividual2.push_back(individual2.at(j));
        }
        
        //for the remaining spots take elements from the other individual that are not already in newIndividual
        for(int j = 0; j < individual2.size(); j++)
        {
            if (find(newIndividual1.begin(), newIndividual1.end(), individual2.at(j)) == newIndividual1.end())
            {
                newIndividual1.push_back(individual2.at(j));
            }
            if (find(newIndividual2.begin(), newIndividual2.end(), individual1.at(j)) == newIndividual2.end())
            {
                newIndividual2.push_back(individual1.at(j));
            }
        }
        
        //get new values
        population.at(pos.at(i)) = newIndividual1;
        population.at(pos.at(i+1)) = newIndividual2;
        
        individual1.clear();
        individual2.clear();
        newIndividual1.clear();
        newIndividual2.clear();
        
    }
}


vector<int> neighbor(vector<int> candidate, int position)
{
    vector<int> neighbor;
    neighbor = candidate;
    
    int point = randomRange(candidate.size());
    while(point == position)
    {
        point = randomRange(candidate.size());
    }
    
    neighbor.at(position) = candidate.at(point);
    neighbor.at(point) = candidate.at(position);
    
    return neighbor;
}

vector<int> Improve(vector<int> solution, vector<vector<float >> coord)
{
    vector<int> candidat;
    vector<float> arg_cand;
    vector<int> best;
    int optiune = 1;

    float min_found = evaluate_solution(coord, solution);
    
    for(int i = 0; i < solution.size() - 1; ++i)
    {
        candidat = neighbor(solution, i);
        
        //calculam functia in candidat
        float fc = evaluate_solution(coord, candidat);
        
        if(fc < min_found && optiune == 0)
        {
            min_found = fc;
            best = candidat;
            break;
        }
        if(fc < min_found && optiune == 1)
        {
            min_found = fc;
            best = candidat;
        }
    }
    if(best.empty())
    {
        return solution;
    }
    //cout<<min_gasit;
    return best;
}

void hill(vector<vector<float >> coord, int restart, vector<int> order, ofstream &file)
{
    
    while(restart)
    {
        vector<int> best;
        unsigned long nrofCities = coord.size();
        if(order.empty()){
            best = generate_solution(nrofCities);
        }
        else
        {
            best = order;
        }
        float best_value = evaluate_solution(coord, best);
        
        int t = 100;
        while(t)
        {
            int local = 1;
            vector<int> candidat ;
            //select candidate solution
            if(order.empty())
            {
                candidat =  generate_solution(nrofCities);
            }
            else
            {
                candidat = order;
            }
            
            
            // evaluate candidate solution
            float candidat_value = evaluate_solution(coord, candidat);

            while(local){
                
                //find a new solution better than candidate
                vector<int> new_sol = Improve(candidat, coord);
                float new_sol_value = evaluate_solution(coord, new_sol);
        
                if(new_sol_value < candidat_value)
                {
                    candidat.clear();
                    candidat = new_sol;
                    candidat_value = new_sol_value;
                }
                else
                {
                    local = 0;
                }
            }
            
            if(candidat_value < best_value)
            {
                best.clear();
                best = candidat;
                best_value = candidat_value;
            }
            t--;
            candidat.clear();
            //cout << best_value << endl;
        }
        file << best_value << endl;
        cout << best_value << "\n";
        restart--;
    }
    
}

void calculateOptimum(string fileName, vector<vector<float >> coord)
{
    vector<int> order;
    ifstream infile(fileName);
    int number;
    while (infile >> number )
    {
        order.push_back(number);
    }
    displaySolution(order);
    for(int i = 0; i < order.size(); i++)
    {
        order.at(i)--;
    }
    displaySolution(order);
    cout << order.size() << endl;
    cout << evaluate_solution(coord, order)<<endl;
    
}

vector<int> TwoOptSwap( const int i, const int k, vector<int> tour)
{
    
    vector<int> new_tour;
    // 1. take route[0] to route[i-1] and add them in order to new_route
    for ( int c = 0; c <= i - 1; ++c )
    {
        new_tour.push_back(tour.at(c));
    }
  
    // 2. take route[i] to route[k] and add them in reverse order to new_route
    int dec = 0;
    for ( int c = i; c <= k; ++c )
    {
        new_tour.push_back(tour.at(k - dec));
        dec++;
    }
  
    // 3. take route[k+1] to end and add them in order to new_route
    for ( int c = k + 1; c < tour.size(); ++c )
    {
        new_tour.push_back(tour.at(c));
    }
    
    return new_tour;
}

// Do all 2-opt combinations
void twoOpt(vector<vector<float> > coord, vector<vector<int> >& population, int posmin)
{

        vector<int> tour = population.at(posmin);
        
        // repeat until no improvement is made
        int improve = 0;
        
        while ( improve < 20 )
        {
            int best_distance = evaluate_solution(coord, tour);
        
            for (int j = 0; j < tour.size() - 1; j++)
            {
                for ( int k = j + 1; k < tour.size(); k++)
                {
                    
                    vector<int> new_tour = TwoOptSwap( j, k, tour);
        
                    int new_distance = evaluate_solution(coord, new_tour);
        
                    if ( new_distance < best_distance )
                    {
                        // Improvement found so reset
                        improve = 0;
                        tour.clear();
                        tour = new_tour;
                        best_distance = new_distance;
                    }
                }
            }
            improve ++;
        }
        population.at(posmin) = tour;
    /*
        //do some random 2optswap for the rest of the population
    for(int i = 0; i < population.size(); i++)
    {
        vector<int> tour = population.at(i);
        
        // repeat until no improvement is made
        int improve = 0;
        
        while ( improve < 20 )
        {
            int best_distance = evaluate_solution(coord, tour);
            
            int nrOfSwaps = 5;
            while(nrOfSwaps)
            {
                int j = randomRange(tour.size() - 2);
                int k = randomInRange(j+1, tour.size() - 1);
                //cout<<j<<" "<<k<< endl;
                if(j < k)
                {
                    vector<int> new_tour = TwoOptSwap( j, k, tour);
                    
                    int new_distance = evaluate_solution(coord, new_tour);
                    
                    if ( new_distance < best_distance )
                    {
                        // Improvement found so reset
                        improve = 0;
                        tour.clear();
                        tour = new_tour;
                        best_distance = new_distance;
                    }
                }
                nrOfSwaps--;
                
            }
            if(nrOfSwaps == 0)
            {
                improve = 20;
            }
        }
        population.at(i) = tour;
    }*/
}

void genetic(vector<vector<float >> coord, int pop_size, ofstream &file)
{
    unsigned long nrofCities = coord.size();
    vector<float> fitness;
    int bestDistance;
    int generations = 1000;
    //generate intial population
    vector<vector<int> > population = generate_population(pop_size, nrofCities);
    
    //calculate fitness
    fitness = calculate_fitness(population, coord);
    
    //get best out of this population
    bestDistance = get_best_value(population, fitness, coord);
    int posmin;
    posmin = get_position_best(fitness);
    vector<int> minFound = population.at(posmin);
    //file << bestDistance << "\n";
    
    while(generations)
    {
        crossover(population, generations);
        mutation(population, generations);
        if(generations <= 1)
        {
            twoOpt(coord, population, posmin);
        }
        //calcultate new fitness
        fitness.clear();
        fitness = calculate_fitness(population, coord);
        
        //selection
        vector<vector<int> > new_population = selection(population, fitness);
        population.clear();
        population = new_population;
        
        //calcultate new fitness
        fitness.clear();
        fitness = calculate_fitness(population, coord);
        
        //evaluate
        float minGeneration = get_best_value(population, fitness, coord);
        cout <<generations << " " << minGeneration << endl;
        //file << minGeneration << "\n";
        
        
        if(bestDistance > minGeneration)
        {
            minFound.clear();
            bestDistance = minGeneration;
            posmin = get_position_best(fitness);
            minFound = population.at(posmin);
        }
        
        generations--;
    }
    file << bestDistance <<"\n";
    cout << bestDistance << endl;
    //displaySolution(minFound);
    //hill(coord, 10, minFound);
}



int main()
{
    clock_t tStart = clock();
    srand(static_cast<unsigned int>(time(NULL)));
    //int pop_size = 100;
    int restart = 20;
    vector<vector<float >> coord;
    vector<int> start;
    //calculateOptimum("att48optimum.txt", coord);
    /*
    for(int i = 0; i < coord.size(); i++)
    {
        cout<<i<< " " <<coord[i].at(0)<<" "<<coord[i].at(1)<<endl;
    }
    
    ofstream f;
    f.open ("berlin52_hill.txt");
    for(int i = 0; i <= 30; i++)
    {
        cout<<i<<" ";
        genetic(coord, pop_size, f);
    }
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    
    f.close();
    
    tStart = clock();
    coord.clear();
    readData("bier127.txt", coord);
    f.open ("bier127_min.txt");
    for(int i = 0; i <= 30; i++)
    {
        cout<<i<<" ";
        genetic(coord, pop_size, f);
    }
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    
    f.close();
    */
    ofstream f;
    f.open ("berlin52_30.txt");
    readData("berlin52.txt", coord);
    for(int i = 0; i <= 30; i++)
    {
        cout<<i<<" ";
        //genetic(coord, pop_size, f);
        hill(coord, restart, start, f);
    }
    //hill(coord, restart, start, f);
    //genetic(coord, pop_size, f);
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    f.close();
    
    tStart = clock();
    f.open ("bier127_30.txt");
    coord.clear();
    readData("bier127.txt", coord);
    for(int i = 0; i <= 30; i++)
    {
        cout<<i<<" ";
        //genetic(coord, pop_size, f);
        hill(coord, restart, start, f);
    }
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    f.close();
    
    tStart = clock();
    f.open ("ch130_30.txt");
    coord.clear();
    readData("ch130.txt", coord);
    for(int i = 0; i <= 30; i++)
    {
        cout<<i<<" ";
        //genetic(coord, pop_size, f);
        hill(coord, restart, start, f);
    }
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    f.close();
    
    tStart = clock();
    f.open ("pr76_30.txt");
    coord.clear();
    readData("pr76.txt", coord);
    for(int i = 0; i <= 30; i++)
    {
        cout<<i<<" ";
        //genetic(coord, pop_size, f);
        hill(coord, restart, start, f);
    }
    f << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
    f.close();
     
    return 0;
}
