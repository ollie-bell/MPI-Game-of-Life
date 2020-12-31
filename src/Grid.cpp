#include "../include/Grid.h"
#include <iostream>
#include <string>
#include <cstring>

Grid::Grid(int n_rows, int n_cols, int id, int iprocs, int jprocs, std::vector<int> all_dims) : rows(n_rows), cols(n_cols), proc_id(id), iprocs(iprocs), jprocs(jprocs), all_dims(all_dims)
{
    curr_gen = nullptr;
    next_gen = nullptr;
}

/*
Initialise grid and allocate memory
*/
void Grid::Initialize()
{
    // allocate memory including padding
    curr_gen = new bool *[rows + 2];
    next_gen = new bool *[rows + 2];
    for (int n = 0; n < rows + 2; n++)
    {
        curr_gen[n] = new bool[cols + 2];
        next_gen[n] = new bool[cols + 2];
    }

    // ensure all cells initialised to false
    for (int i = 0; i < rows + 2; i++)
        for (int j = 0; j < cols + 2; j++)
        {
            curr_gen[i][j] = false;
            next_gen[i][j] = false;
        }

    // initialise MPI I/O
    init_MPI_IO();
}

/*
Fill the grid with a random state given a probability [0, 1] of alive cells
*/
void Grid::Randomize(double prob_alive)
{
    Initialize();

    // padding excluded from loop
    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++)
        {
            double val = (double)((double)rand() / (double)RAND_MAX);
            if (val <= prob_alive)
                curr_gen[i][j] = true;
        }
}

/*
Load initial state from binary input files
*/
void Grid::from_infile(const char *fname)
{
    Initialize();
    char *temp = new char[rows * cols * sizeof(bool)];

    MPI_File fh;
    MPI_Offset disp = 0;
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, disp, MPI_BYTE, filetype, "native", MPI_INFO_NULL);
    MPI_File_read_all(fh, temp, 1, datatype, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    for (int i = 1; i <= rows; i++)
        memcpy(&curr_gen[i][1], &temp[(i - 1) * cols], sizeof(bool) * cols);

    delete[] temp;
}

/*
Copy the grid state into a buffer and update the size of data in the buffer
*/
void Grid::write_to_buffer(char *buffer, int grid_size, int &buff_size)
{
    for (int i = 1; i <= rows; i++)
        memcpy(&buffer[buff_size + (i - 1) * cols], &curr_gen[i][1], sizeof(bool) * cols);
    buff_size += grid_size;
}

/*
Dump the buffer into a binary file and reset the buffer
*/
void Grid::dump_buffer(char *buffer, int &buff_size, int &dump_cnt, int freq_dump)
{
    int num_dump = buff_size / (rows * cols * sizeof(bool));

    if (proc_id == 0)
        std::cout << "Dumping grid states   " << dump_cnt * freq_dump << " - "
                  << dump_cnt * freq_dump + buff_size / (rows * cols * sizeof(bool)) - 1
                  << "   to file" << std::endl;

    std::string fname = "./outfiles/dump" + std::to_string(dump_cnt);

    MPI_File fh;
    MPI_Offset disp = 0;
    MPI_File_open(MPI_COMM_WORLD, fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view(fh, disp, MPI_BYTE, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, buffer, num_dump, datatype, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    buff_size = 0;
    dump_cnt++;
}

/*
Initialise filetype and datatype for MPI I/O
*/
void Grid::init_MPI_IO()
{
    int start_i = 0;
    int i_id = proc_id / jprocs;
    for (int i = 0; i < i_id; i++)
        start_i += all_dims[2 * i];

    int start_j = 0;
    int j_id = proc_id % jprocs;
    for (int j = 0; j < j_id; j++)
        start_j += all_dims[2 * j + 1];

    int n = all_dims.size();
    const int gsize_array[2] = {all_dims[n - 2], all_dims[n - 1]};
    const int gsubsize_array[2] = {rows, cols};
    const int gstart_array[2] = {start_i, start_j};

    MPI_Type_create_subarray(2, gsize_array, gsubsize_array, gstart_array, MPI_ORDER_C, MPI_BYTE, &filetype);
    MPI_Type_commit(&filetype);

    MPI_Type_contiguous(rows * cols, MPI_BYTE, &datatype);
    MPI_Type_commit(&datatype);
}

/*
Execute all communications with neighbouring processors to populate padding cells
*/
void Grid::do_comms(MPI_Request *requests, int &req_cnt, int n)
{
    // alternates between custom MPI_Datatype's for the curr_gen and next_gen array's
    req_cnt = 0;
    for (size_t k = 0; k < all_comms.size() / 2; k++)
    {
        int send_ind = 2 * k + n % 2;
        int recv_ind = (all_comms.size() - 2) - 2 * k + n % 2; // the recvs are stored in the reverse order relative to the sends

        MPI_Irecv(MPI_BOTTOM, 1, all_comms[recv_ind].comm_types[1], all_comms[recv_ind].neighbour_proc, 0, MPI_COMM_WORLD, &requests[req_cnt]);
        req_cnt++;
        MPI_Isend(MPI_BOTTOM, 1, all_comms[send_ind].comm_types[0], all_comms[send_ind].neighbour_proc, 0, MPI_COMM_WORLD, &requests[req_cnt]);
        req_cnt++;
    }
}

/*
Determine whether a cell will be alive or dead in the next generation
*/
void Grid::next_gen_state(int i, int j, bool curr_state)
{
    // count number of alive neighbours
    int count = 0;
    for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
        {
            int ni = i + ii;
            int nj = j + jj;
            count += curr_gen[ni][nj]; // using implicit bool to int conversion
        }
    count -= curr_state; // don't count itself

    // apply rules based on current state and number of alive neighbours
    if (!curr_state && count != 3)
        next_gen[i][j] = false;
    else if (count == 3)
        next_gen[i][j] = true;
    else if (curr_state && count != 2)
        next_gen[i][j] = false;
    else
        next_gen[i][j] = true;
}

/*
Swap the current generation grid and the next generation grid
*/
void Grid::update_grid()
{
    bool **tmp = curr_gen;
    curr_gen = next_gen;
    next_gen = tmp;
}

/*
Play Conway's Game of Life on the grid with or without periodic boundaries. 

The grid state is copied to a buffer after each generation, and the buffer is 
dumped to a binary file every freq_dump generations.
*/
void Grid::Game_of_Life(int num_gen, int freq_dump)
{
    // initialising writing to buffer and output file dump
    int dump_cnt = 0;
    int buff_size = 0;
    int grid_size = rows * cols * sizeof(bool);
    int max_buff_size = grid_size * freq_dump;
    char *buffer = new char[max_buff_size];
    // init_MPI_IO();

    // write initial state to buffer
    write_to_buffer(buffer, grid_size, buff_size);
    if (buff_size == max_buff_size)
        dump_buffer(buffer, buff_size, dump_cnt, freq_dump);

    // set-up requests for MPI communications
    int req_cnt;
    MPI_Request *requests = new MPI_Request[all_comms.size()];

    // *** SIMULATION START *** //
    for (int n = 0; n < num_gen; n++)
    {
        // execute all communications at boundaries between sub-grids
        do_comms(requests, req_cnt, n);
        MPI_Waitall(req_cnt, requests, MPI_STATUSES_IGNORE); // waiting outside do_comms for flexibility on when the Waitall is called
                                                             // e.g. could update cells (2:rows-1, 2:cols-1) before Waitall, then finish updating the edge cells

        // To each cell apply the game rules and determine its state for the next generation
        for (int i = 1; i <= rows; i++)
            for (int j = 1; j <= cols; j++)
                next_gen_state(i, j, curr_gen[i][j]);

        // Swap the current and next generation grids
        update_grid();

        // output grid data
        write_to_buffer(buffer, grid_size, buff_size);
        if (buff_size == max_buff_size)
        {
            dump_buffer(buffer, buff_size, dump_cnt, freq_dump);
        }
    }
    // *** SIMULATION FINISH *** //

    // final file dump on remaining buffer items
    if (buff_size > 0)
        dump_buffer(buffer, buff_size, dump_cnt, freq_dump);

    // tidy up before exiting
    delete[] buffer;
    delete[] requests;

    MPI_Type_free(&filetype);
    MPI_Type_free(&datatype);

    for (size_t i = 0; i < all_comms.size(); i++)
        all_comms[i].free_types();
}
