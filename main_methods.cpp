// C1
int main()
{
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		
        NearestNeighbour(0);
		UpdateMinTour();

		OutputMinTour();
    }

	return 0;
}

//C2
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
	
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		
		auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            NearestNeighbour(rand_idx);
            UpdateMinTour();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

		OutputMinTour();
    }

    return 0;
}

// C3
int main()
{
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
		
        Greedy();
		UpdateMinTour();

		OutputMinTour();
    }

	return 0;
}

//C4
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
	
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
		computeClosestNeighbours();
		
        Greedy();
		UpdateMinTour();
		
		auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            GreedyRandomNode(rand_idx);
            UpdateMinTour();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

		OutputMinTour();
    }

    return 0;
}

//C5
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
	
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
		computeClosestNeighbours();
		
        Greedy();
		UpdateMinTour();
		
		auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            GreedyRandomShuffle();
            UpdateMinTour();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

		OutputMinTour();
    }

    return 0;
}

// C6
int main()
{
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		addAllNodes();
		
        GreedyTour(1);
		UpdateMinTour();

		OutputMinTour();
    }

	return 0;
}

//C7
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
	
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		addAllNodes();
		
		auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            GreedyTourRandom();
            UpdateMinTour();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

		OutputMinTour();
    }

    return 0;
}

//C8
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
	
	readInput();
	
	if (n <= 3) // handle trivial corner cases
    {
		for (int i = 0; i < n; i++)
			cout << i << '\n';
    }
    else {
		computeDistances(); 
		num_edges = (n*(n - 1)) / 2;
		
		auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

			int rand_idx = (rand() % (n - 1)) + 1;
            ClarkeWright(rand_idx);
            UpdateMinTour();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

		OutputMinTour();
    }

    return 0;
}

//C9
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
        computeClosestNeighboursTwoOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            NearestNeighbour(rand_idx);
            TwoOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C10
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		addAllNodes();
        computeClosestNeighboursTwoOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            GreedyTourRandom();
            TwoOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C11
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
		computeClosestNeighbours();
        computeClosestNeighboursTwoOpt();
        
        Greedy();
        TwoOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

			int rand_idx = (rand() % (n - 1)) + 1;
            GreedyRandomNode(rand_idx);
            TwoOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C12
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
        num_edges = (n*(n - 1)) / 2;
        sortEdges();
        computeClosestNeighboursTwoOpt();
        
        Greedy();
        TwoOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            GreedyRandomShuffle();
            TwoOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C13
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
        computeClosestNeighboursTwoOpt();
        
        Greedy();
        TwoOpt();

        OutputMinTour();
    }

    return 0;
}

//C14
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		num_edges = (n*(n - 1)) / 2;
        computeClosestNeighboursTwoOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();
			
			int rand_idx = (rand() % (n - 1)) + 1;
            ClarkeWright(rand_idx);
            TwoOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C15
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
        computeClosestNeighbours();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1980000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            NearestNeighbour(rand_idx);
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C16
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		addAllNodes();
        computeClosestNeighbours();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1980000) {
            auto start_run = high_resolution_clock::now();

            GreedyTourRandom();
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C17
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
		computeClosestNeighbours();
        
        Greedy();
        ThreeOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

			int rand_idx = (rand() % (n - 1)) + 1;
            GreedyRandomNode(rand_idx);
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C18
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
        num_edges = (n*(n - 1)) / 2;
        sortEdges();
        computeClosestNeighbours();
        
        Greedy();
        ThreeOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            GreedyRandomShuffle();
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C19
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		num_edges = (n*(n - 1)) / 2;
		sortEdges();
        computeClosestNeighbours();
        
        Greedy();
        ThreeOpt();

        OutputMinTour();
    }

    return 0;
}

//C20
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
		num_edges = (n*(n - 1)) / 2;
        computeClosestNeighbours();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();
			
			int rand_idx = (rand() % (n - 1)) + 1;
            ClarkeWright(rand_idx);
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C21
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        computeDistances();
        computeClosestNeighboursQuadrants();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1980000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            NearestNeighbour(rand_idx);
            ThreeOptQuadrants();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C22
int main()
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
        
    computeDistances();
    
    if (n <= 3) // handle trivial corner cases
    {
        for (int i = 0; i < n; i++)
            cout << i << '\n';
    }
    else {
        num_edges = (n*(n - 1)) / 2;
        sortEdges();
        addAllNodes();
        computeClosestNeighbours();
        
        GreedyTourRandom();
        ThreeOpt();
        Greedy();
        ThreeOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            NearestNeighbour(rand_idx);
            ThreeOpt();
            
            rand_idx = (rand() % (n - 1)) + 1;
            ClarkeWright(rand_idx);
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C23
{
    high_resolution_clock::time_point time = high_resolution_clock::now();
    
    readInput();
        
    computeDistances();
    
    if (n <= 17) // handle trivial corner cases
    {
        DynProg();
        
        OutputMinTour();
    }
    else {
        num_edges = (n*(n - 1)) / 2;
        sortEdges();
        addAllNodes();
        computeClosestNeighbours();
        
        GreedyTourRandom();
        ThreeOpt();
        Greedy();
        ThreeOpt();
        
        auto duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count(); // duration will be time passed
        auto last_run = duration; // last_run will be time last iteration took
        srand(std::time(0));
        while (duration + (1.05 * last_run) < 1990000) {
            auto start_run = high_resolution_clock::now();

            int rand_idx = (rand() % (n - 1)) + 1;
            NearestNeighbour(rand_idx);
            ThreeOpt();
            
            rand_idx = (rand() % (n - 1)) + 1;
            ClarkeWright(rand_idx);
            ThreeOpt();

            last_run = duration_cast<microseconds>(high_resolution_clock::now() - start_run).count();
            duration = duration_cast<microseconds>(high_resolution_clock::now() - time).count();
        }

        OutputMinTour();
    }

    return 0;
}

//C24
int main()
{
	readInput();

	if (n <= 3) // handle trivial corner cases
	{
		for (int i = 0; i < n; i++)
			cout << i << '\n';
	}
	else {
		computeDistances();
		num_edges = (n*(n - 1)) / 2;

		ClarkeWright(0);
		UpdateMinTour();

		OutputMinTour();
	}

	return 0;
}