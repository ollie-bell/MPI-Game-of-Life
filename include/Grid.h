#ifndef GRID_H
#define GRID_H
#include <mpi.h>
#include <vector>
#include "Comms.h"

/*
Class to store game of life grid data and methods for playing the game of life
*/
class Grid
{
public:
    // constructor & destructor
    Grid(int n_rows, int n_cols, int id, int iprocs, int jprocs, std::vector<int> all_dims);
    ~Grid()
    {
        for (int n = 0; n < rows + 2; n++)
            delete[] curr_gen[n];
        delete[] curr_gen;

        for (int n = 0; n < rows + 2; n++)
            delete[] next_gen[n];
        delete[] next_gen;
    };

    // Grid variables
    int rows, cols; // exlcusive of padding
    int proc_id;
    bool **curr_gen;
    bool **next_gen;

    // Global grid information
    int iprocs;
    int jprocs;
    std::vector<int> all_dims;

    // container for all communications with neighbouring processors
    std::vector<Comms> all_comms;

    // Allocate memory and initialise MPI-I/O
    void Initialize();
    // Fill the grid with a random state given a probability [0, 1] of alive cells
    void Randomize(double prob_alive);
    // Load initial state from binary input files
    void from_infile(const char *fname);

    // Play Conway's Game of Life on the grid with or without periodic boundaries
    void Game_of_Life(int num_gen, int freq_dump);

private:
    MPI_Datatype filetype;
    MPI_Datatype datatype;

    // Copy the grid state into a buffer and update the size of data in the buffer
    void write_to_buffer(char *buffer, int grid_size, int &buff_size);
    // Dump the buffer into a binary file and reset the buffer
    void dump_buffer(char *buffer, int &buff_size, int &dump_cnt, int freq_dump);
    // Initialise filetype and datatype for MPI I/O
    void init_MPI_IO();

    // Execute all communications with neighbouring processors to populate padding cells
    void do_comms(MPI_Request *requests, int &req_cnt, int n);
    // Determine whether a cell will be alive or dead in the next generation
    void next_gen_state(int i, int j, bool curr_state);
    // Swap the current generation grid and the next generation grid
    void update_grid();
};

#endif