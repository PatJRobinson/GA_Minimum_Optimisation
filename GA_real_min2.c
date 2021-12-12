/*
This program is for the 1st function defined in the assignment brief. 
It tries to find the solution that yeilds the minimum value, using a genetic algorithm.
Function = sum((x[i+1] - x^2)^2 + (1 - x)^2), 
for values 0 <= x < 20, where x is 
a gene ranging from -100 to 100.
*/


#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define P 50
#define N 20
#define MAX 5.12
#define MIN -5.12
#define ITERATIONS 50
#define FUNCTION_CHOICE 1 // 0, 1, 2

float MUTRATE = 0.05;
float MUT_ALTER = 50;

FILE * solution_fpt;
FILE * results_fpt;

typedef struct{
    float gene[N];
    float fitness;
} individual;

// keep track of the best solution so far
individual current_best;

// a function which tests an individuals fitness and returns it as an integer value
float test_function(individual ind){

    float utility = 0;

    
    for (int i = 0; i < N-1; i++){

        float gene_fitness;
        if (FUNCTION_CHOICE == 0) 
        {
            //sum((x[i+1] - x^2)^2 + (1 - x)^2)
            gene_fitness = 100 * (pow(ind.gene[i+1] - pow(ind.gene[i], 2), 2) + pow(1-ind.gene[i], 2));
            utility += gene_fitness;
        }
        else if (FUNCTION_CHOICE == 1)
        {
            /*10*N + sum(x^2 - 10.cos(2*PI*x))----------------------*/
            gene_fitness = pow(ind.gene[i], 2) - 10 * cos(2*M_PI*ind.gene[i]);
            utility += gene_fitness;
        }

    }
    if (FUNCTION_CHOICE == 1)
       utility += 10*N;
    
    return utility;
}

individual population[P];
// offspring population
individual offspring[P];

void seed(void) {
    //1. randomly seed a population of individuals

    for (int i = 0; i< P; i++){
        for (int j = 0; j < N; j++){
            // setting each gene to a random binary digit
            if (FUNCTION_CHOICE == 0)
                population[i].gene[j] = rand()%(int)((MAX - MIN)* 100)/100.0 + MIN;//roundf(((rand()%20001/100.0) - 100) * 100) / 100;//(rand()%199) - 100 + (rand()%1001/1000.0);
            else if (FUNCTION_CHOICE == 1)
                population[i].gene[j] = rand()%(int)((MAX - MIN)* 100)/100.0 + MIN; // -5.12 < x < 5.12
        }
        population[i].fitness = test_function(population[i]);
    }
}


// a function which applies tournament selection, crossover and mutation to produce a new generation of offspring 
void iterate(individual *curr_best, int iteration) {
    
    //2. tournament selection - improvement over roulette wheel which can 'prematurely converge' onto local optima
    for (int i = 0; i < P; i++){
        int parent1 = rand()%P;
        int parent2 = rand()%P;
        if (population[parent1].fitness < population[parent2].fitness)
            offspring[i] = population[parent1];
        else   
            offspring[i] = population[parent2];
    }

    // initialise best ever solution on first iteration
    if (iteration == 0)
        current_best.fitness = test_function(offspring[0]);
    // Need to keep a copy of the best ever solution
    int worst_at = 0;
    individual worst = offspring[0];
    individual best = offspring[0];
    best.fitness = test_function(offspring[0]);
    for (int i = 1; i < P; i++) {

        offspring[i].fitness = test_function(offspring[i]);
        if (offspring[i].fitness < best.fitness){
            best = offspring[i];
        }
        if (offspring[i].fitness > worst.fitness){
            worst = offspring[i];
            worst_at = i;
        }
    }
    if (best.fitness > current_best.fitness){
        offspring[worst_at] = current_best;
    } else {
        current_best = best;
    }
    
    //3. Crossover

    individual temp;

    for (int i = 0; i<P; i+=2){

        temp = offspring[i];
        // swapping tails
        int crosspoint = rand() % N;
        for (int j = crosspoint; j<N ; j++){
            offspring[i].gene[j] = offspring[i+1].gene[j];
            offspring[i+1].gene[j] = temp.gene[j];
        }

    }

    //3. Mutation
    for (int i = 0; i < P; i++){
        for (int j = 0; j < N; j++){
            // mutation occurs
            float mut_prob = rand()%1001/1000.0;
            if (mut_prob < MUTRATE) {
                if (mut_prob < MUTRATE/2) offspring[i].gene[j] += MUT_ALTER;
                else offspring[i].gene[j] -= MUT_ALTER;

                if (offspring[i].gene[j] > 1) offspring[i].gene[j] = 1;
                else if (offspring[i].gene[j] < 0) offspring[i].gene[j] = 0;
                
            }
        }
        //4. re-test the offspring population
        offspring[i].fitness = test_function(offspring[i]);
    }

    
    
}

/*
A function to print solution data to a csv file
*/
void print_all_solutions_to_csv(int iteration, char parentOrOffspring){

    individual * pop_ptr;
    if (parentOrOffspring == 'P')
        pop_ptr = &population[0];
    else pop_ptr = &offspring[0];

    fprintf(solution_fpt, "Population: %d\n", iteration);
    fprintf(solution_fpt,"\n");
    fprintf(solution_fpt, "Curr Best: \n");
    for (int y = 0; y < N; y++)
        fprintf(solution_fpt, "%f ,", current_best.gene[y]);
    fprintf(solution_fpt, "\n"); 
    fprintf(solution_fpt, "\n");       

    for (int i = 0; i < P; i++){
        for (int y = 0; y < N; y++)
            fprintf(solution_fpt, "%f ,", (pop_ptr+i)->gene[y]);
        fprintf(solution_fpt, "\n");
    }
    fprintf(solution_fpt, "\n");
}

/*
A function to print results data to a csv file
*/
void print_results_to_csv(int iteration, float ave){

    fprintf(results_fpt, "%d: Best fitness = %f     Average fitness = %f\n", iteration, current_best.fitness, ave);

}

int main(void) {


    // seeding the rng
    srand( time(NULL));
    
    seed();

    float mean_offspring_fitness[ITERATIONS];
    float best_offspring_fitness[ITERATIONS];
    // initialise the best solution so far

    /*
    initialize csv files
    */
    // a file to store the all solutions of all populations
    solution_fpt = fopen("Solutions.csv", "w+");
    fprintf(solution_fpt, "All solutions from each population\n");

    // a file to store the mean and best fitnesses of all populations
    results_fpt = fopen("Results.csv", "w+");
    fprintf(results_fpt, "Best and average from each population\n");

    //a debug file
    FILE * debug = fopen("Debug.csv", "w+");
    fprintf(debug, "Best and average from each population\n");

    //print the first parent population
    print_all_solutions_to_csv(0, 'P');
    print_results_to_csv(0, 0);

    for (int x = 0; x < ITERATIONS; x++) {
    
        /* 
        apply tournament selection on parent population, then crossover 
        and mutation on the offspring. Retains best ever solution
        */
        iterate(&current_best, x);

        for (int i = 0; i< P; i++){
            fprintf(debug, "%f", population[i].fitness);

        fprintf(debug, "\n");
        fprintf(debug, "\n");
        }

        /*
        copy offspring population accross to parent population
        and count total fitness for producing average
        */
        int total_offspring_fitness = 0;
        for (int i = 0; i < P; i++) {
            total_offspring_fitness += offspring[i].fitness;
            population[i] = offspring[i];
        }
        float ave = total_offspring_fitness/P;

        printf("Best fitness = %f, mean fitness = %f \n ", current_best.fitness, ave);
        //print offspring population
        print_all_solutions_to_csv(x + 1, 'O');
        print_results_to_csv(x + 1, ave);

    }
}
