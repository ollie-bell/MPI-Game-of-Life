#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <mpi.h>

#include "../include/Grid.h"
#include "../include/Comms.h"

int id, p;

/*
Decompose the global grid across the available processes into as close to square as possible sub-grids
*/
void domain_decomposition(int rows, int cols, int &imax, int &jmax, int &proc_rows, int &proc_cols, std::vector<int> &all_dims)
{
    // loop to determine the best arrangement of processes by checking which pair of factors
    // of p results in sub-grids closest to square
    int gap = rows;
    if (cols > gap)
        gap = cols;

    for (int n = 1; n <= p; n++)
        if (p % n == 0)
        {
            int tmp = abs(rows / n - cols / (p / n));
            if (tmp < gap)
            {
                imax = n;
                jmax = p / n;
                gap = tmp;
            }
        }

    // (i, j) position of this process
    int j_id = id % jmax;
    int i_id = id / jmax;

    // calculate the number of rows in the sub-grid this process is repsonsible for
    int rows_rem = rows;
    for (int i = 0; i <= i_id; i++)
    {
        proc_rows = rows_rem / (imax - i);
        rows_rem -= proc_rows;
    }

    // calculate the number of cols in the sub-grid this process is repsonsible for
    int cols_rem = cols;
    for (int j = 0; j <= j_id; j++)
    {
        proc_cols = cols_rem / (jmax - j);
        cols_rem -= proc_cols;
    }

    // share grid dimensions between processes
    std::vector<int> dims{proc_rows, proc_cols};
    MPI_Allgather(dims.data(), 2, MPI_INT, all_dims.data(), 2, MPI_INT, MPI_COMM_WORLD);
    all_dims[2 * p] = rows;
    all_dims[2 * p + 1] = cols;
}

/*
Set-up all the communication objects necessary to perform sends/receives to/from neighbouring processes
*/
void create_communications(Grid &sub_grid, int imax, int jmax, bool is_periodic)
{
    // (i, j) position of this process
    int j_id = id % jmax;
    int i_id = id / jmax;

    // loop over all potential neighbour processes
    for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
        {
            // (i, j) position of the potential neighbour process
            int i, j;
            if (is_periodic)
            {
                i = (ii + i_id + imax) % imax;
                j = (jj + j_id + jmax) % jmax;
            }
            else
            {
                i = (ii + i_id);
                j = (jj + j_id);
            }

            // ensure assumed neighbour is valid (i.e. not out of bounds, not self)
            if ((i >= 0 && i < imax) && (j >= 0 && j < jmax) && (ii != 0 || jj != 0))
            {
                // process id of the neighbour process
                int n_id = i * jmax + j;

                // start and end indices of grid data to be sent and received
                int start_i_send, start_j_send, end_i_send, end_j_send;
                int start_i_recv, start_j_recv, end_i_recv, end_j_recv;

                if (ii == -1 && jj == -1) // --> top left corner send/recv
                {
                    start_i_send = 1, start_j_send = 1;
                    end_i_send = 1, end_j_send = 1;
                    start_i_recv = 0, start_j_recv = 0;
                    end_i_recv = 0, end_j_recv = 0;
                }
                else if (ii == -1 && jj == 0) // --> top row send/recv
                {
                    start_i_send = 1, start_j_send = 1;
                    end_i_send = 1, end_j_send = sub_grid.cols;
                    start_i_recv = 0, start_j_recv = 1;
                    end_i_recv = 0, end_j_recv = sub_grid.cols;
                }
                else if (ii == -1 && jj == 1) // --> top right corner send/recv
                {
                    start_i_send = 1, start_j_send = sub_grid.cols;
                    end_i_send = 1, end_j_send = sub_grid.cols;
                    start_i_recv = 0, start_j_recv = sub_grid.cols + 1;
                    end_i_recv = 0, end_j_recv = sub_grid.cols + 1;
                }
                else if (ii == 0 && jj == -1) // --> left column send/recv
                {
                    start_i_send = 1, start_j_send = 1;
                    end_i_send = sub_grid.rows, end_j_send = 1;
                    start_i_recv = 1, start_j_recv = 0;
                    end_i_recv = sub_grid.rows, end_j_recv = 0;
                }
                else if (ii == 0 && jj == 1) // --> right column send/recv
                {
                    start_i_send = 1, start_j_send = sub_grid.cols;
                    end_i_send = sub_grid.rows, end_j_send = sub_grid.cols;
                    start_i_recv = 1, start_j_recv = sub_grid.cols + 1;
                    end_i_recv = sub_grid.rows, end_j_recv = sub_grid.cols + 1;
                }
                else if (ii == 1 && jj == -1) // --> bottom left corner send/recv
                {
                    start_i_send = sub_grid.rows, start_j_send = 1;
                    end_i_send = sub_grid.rows, end_j_send = 1;
                    start_i_recv = sub_grid.rows + 1, start_j_recv = 0;
                    end_i_recv = sub_grid.rows + 1, end_j_recv = 0;
                }
                else if (ii == 1 && jj == 0) // --> bottom row send/recv
                {
                    start_i_send = sub_grid.rows, start_j_send = 1;
                    end_i_send = sub_grid.rows, end_j_send = sub_grid.cols;
                    start_i_recv = sub_grid.rows + 1, start_j_recv = 1;
                    end_i_recv = sub_grid.rows + 1, end_j_recv = sub_grid.cols;
                }
                else if (ii == 1 && jj == 1) // --> bottom right corner send/recv
                {
                    start_i_send = sub_grid.rows, start_j_send = sub_grid.cols;
                    end_i_send = sub_grid.rows, end_j_send = sub_grid.cols;
                    start_i_recv = sub_grid.rows + 1, start_j_recv = sub_grid.cols + 1;
                    end_i_recv = sub_grid.rows + 1, end_j_recv = sub_grid.cols + 1;
                }

                // Define a temporary communications object for communications with the current neighbour
                Comms temp;
                temp.neighbour_proc = n_id;
                // store the start and end indices for where grid data to be sent is and where to receive
                // grid data into for the this process
                temp.comm_inds[0][0] = start_i_send;
                temp.comm_inds[0][1] = start_j_send;
                temp.comm_inds[0][2] = end_i_send;
                temp.comm_inds[0][3] = end_j_send;
                temp.comm_inds[1][0] = start_i_recv;
                temp.comm_inds[1][1] = start_j_recv;
                temp.comm_inds[1][2] = end_i_recv;
                temp.comm_inds[1][3] = end_j_recv;
                // create MPI types for these communications
                // and store the Comms object in a vector held in the Grid object for this process
                temp.create_MPI_types(sub_grid.curr_gen);
                sub_grid.all_comms.push_back(temp);
                temp.create_MPI_types(sub_grid.next_gen); // swapping between curr_gen and next_gen, which have different memory locations
                sub_grid.all_comms.push_back(temp);       // therefore need to define a seperate MPI type for communication using these arrays
            }
        }
}

/*
Write out useful information about the simulation for automated post-processing
*/
void write_metadata(int rows, int cols, int iprocs, int jprocs, bool is_periodic, int num_gen, int freq_dump, std::vector<int> all_dims)
{
    if (id == 0)
    {
        std::ofstream metadata;
        metadata.open("./outfiles/metadata.txt", std::ios::out);

        metadata << "grid_size: " << rows << " " << cols << std::endl;
        metadata << "num_procs: " << p << std::endl;
        metadata << "proc_layout: " << iprocs << " " << jprocs << std::endl;
        metadata << "num_gen: " << num_gen << std::endl;
        metadata << "periodic: " << is_periodic << std::endl;
        metadata << "dump_freq: " << freq_dump << std::endl;

        for (int k = 0; k < p; k++)
            metadata << "sub_grid_" << k << ": " << all_dims[2 * k] << " " << all_dims[2 * k + 1] << std::endl;

        metadata.close();
    }
}

void show_usage(const char *program_name)
{
    if (id == 0)
    {
        std::cerr << "Usage:\tmpiexec [--use-hw-threads] -np <#> " << program_name << " <args> [opts]\n"
                  << std::endl
                  << "args :"
                  << "\tnrows (unsigned int) --> number of rows in the grid\n"
                  << "\tncols (unsigned int) --> number of columns in the grid\n"
                  << "\tngens (unsigned int) --> number of generations to simulate\n"
                  << "\tEITHER filename (string) OR prob_alive (float) --> initial state of the grid\n"
                  << std::endl
                  << "opts :"
                  << "\t--periodic --> enable periodic boundaries [default: disabled]\n"
                  << "\t--dump=<unsigned int> --> number of generations between output file dumps [default: 101]\n"
                  << std::endl
                  << "e.g. :"
                  << "\tmpiexec -np 4 " << program_name << " 247 331 250 grid.bin --periodic --dump=63\n"
                  << std::endl
                  << "Please refer to the README for further information\n"
                  << std::endl;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    srand(time(NULL) + id * 42);

    // input definitions
    int rows, cols, num_gen, freq_dump;
    bool is_periodic, from_file;
    double prob_alive;
    const char *filename;

    // parse command line arguments
    if (argc < 5 || argc > 7)
    {
        show_usage(argv[0]);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 1; // FIXME: MPI prints to stdout informing about non-zero exit code - how do I exit gracefully?
    }
    else
    {
        // determine if user has input a filename or probability for the initial state
        if (id == 0)
        {
            char *endptr;
            double input_val = strtod(argv[4], &endptr);

            if (*endptr != '\0')
            {
                std::ifstream file_exists(argv[4]);
                if (file_exists.fail())
                {
                    std::cerr << "File " << argv[4] << " not found\n"
                              << std::endl;
                    file_exists.close();
                    MPI_Abort(MPI_COMM_WORLD, 1); // FIXME: need to exit gracefully without using MPI_Abort
                }
                else
                {
                    from_file = true;
                    file_exists.close();
                }
            }
            else
            {
                // check the value of input_val is between 0.0 and 1.0
                if (input_val < 0.0 || input_val > 1.0)
                {
                    std::cerr << "Initial state specified by a probability of alive cells must be between 0.0 and 1.0\n"
                              << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1); // FIXME: need to exit gracefully without using MPI_Abort
                }
                from_file = false;
            }
        }
        MPI_Bcast(&from_file, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD); // this is blocking

        // args
        rows = atoi(argv[1]);
        cols = atoi(argv[2]);
        num_gen = atoi(argv[3]);
        if (from_file)
        {
            filename = argv[4];
        }
        else
        {
            prob_alive = atof(argv[4]);
        }

        // opts -- default values
        is_periodic = false;
        freq_dump = 101;

        // opts -- override if present
        int i = 5;
        while (i < argc)
        {
            std::string opt(argv[i]);

            if (opt == "--periodic")
            {
                is_periodic = true;
            }
            else if (opt.substr(0, 6) == "--dump")
            {
                freq_dump = std::stoi(opt.substr(7));
            }
            else
            {
                show_usage(argv[0]);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Finalize();
                return 1; // FIXME: MPI prints to stdout informing about non-zero exit code - how do I exit gracefully?
            }
            i++;
        }
    }

    // work out decomposition of the grid over available processes
    int iprocs, jprocs, proc_rows, proc_cols;
    std::vector<int> all_dims((p + 1) * 2);
    domain_decomposition(rows, cols, iprocs, jprocs, proc_rows, proc_cols, all_dims);

    // Create a Grid object for each process
    Grid sub_grid(proc_rows, proc_cols, id, iprocs, jprocs, all_dims);

    // Initialise the grid
    if (from_file)
    {
        sub_grid.from_infile(filename);
    }
    else
    {
        sub_grid.Randomize(prob_alive);
    }

    // Set up communications between processes
    create_communications(sub_grid, iprocs, jprocs, is_periodic);

    if (id == 0)
    {
        std::cout << "-----------------------------" << std::endl;
        std::cout << "  MPI Conway's Game of Life  " << std::endl;
        std::cout << "-----------------------------" << std::endl;
        std::cout << std::endl;

        std::cout << "Grid size: " << rows << " x " << cols << std::endl;
        std::cout << "Process layout: " << iprocs << " x " << jprocs << std::endl;
        std::cout << "Periodic boundaries: " << std::boolalpha << is_periodic << std::endl;

        std::cout << std::endl;

        std::cout << "Simulation running ... " << std::endl;
    }
    // Write out simulation info for post-processing
    write_metadata(rows, cols, iprocs, jprocs, is_periodic, num_gen, freq_dump, all_dims);

    // Run the simulation
    sub_grid.Game_of_Life(num_gen, freq_dump);

    if (id == 0)
        std::cout << "Simulation complete.\n"
                  << std::endl;

    MPI_Finalize();
    return 0;
}
