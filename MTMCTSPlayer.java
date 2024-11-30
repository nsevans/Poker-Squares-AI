import java.util.Random;

import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * A Multi-threaded Monte Carlo Tree Search (MTMCTS) Algorithm with a depth limit. Each thread is used to calculate a score for each
 * possible position a card can be played at in any given turn. The threads all contain their own total and average scores
 * which are used for calculating the best position to play on a given turn.
 * 
 * Author: Nicholas Evans
 */
public class MTMCTSPlayer implements PokerSquaresPlayer{
    
	private static final int SIZE = PokerSquares.SIZE;			// Height and width of the card grid
	private static final int NUM_POS = SIZE * SIZE;				// Number of possible positions in the car grid
	private static final int DEFAULT_DEPTH_LIMIT = 14;			// Default max depth the algorithm can go to
	private static final boolean DEFAULT_VERBOSITY = false;		// Default setting for verbose printing
	private static final int DEFAULT_THREAD_MULTIPLIER = 2;		// Default setting for number of threads run per turn (available positions multiplied by this)

	private Random random;
	private PokerSquaresPointSystem system;
	
	private ArrayList<Card> cardsInDeck;			// Array of all cards in the deck, used to determine what cards are left
	private ArrayList<int[]> availablePositions;	// 2D array to hold the coords of each available position in the grid
	private Card[][] grid;							// Playing grid for the player
    
	private int currentRound;		// Current round of the game
	private long totalPlayTime;		// Tracks the total running time of the algorithm

	private int depthLimit;			// Depth limit for the expectimax function
	private int threadMultiplier;	// Setting for number of threads run per turn (available positions multiplied by this)
	private boolean verbose;		// Setting for verbose printing, used for testing

	/**
	 * Creates a Multi Threaded Monte Carlo Tree Search player that simulates random plays on separate threads
	 * and estimates the best position on the game grid with a default depth limit and a default verbosity.
	 */
	public MTMCTSPlayer()
	{
		this(DEFAULT_DEPTH_LIMIT, DEFAULT_THREAD_MULTIPLIER, DEFAULT_VERBOSITY);
	}

	/**
	 * Creates a MTMCTS player that simulates random plays on separate threads
	 * and estimates the best position on the game grid with a given depth limit and a default verbosity.
	 * @param depthLimit The max depth limit for each simulated play
	 */
	public MTMCTSPlayer(int depthLimit)
	{
		this(depthLimit, DEFAULT_THREAD_MULTIPLIER, DEFAULT_VERBOSITY);
	}

	/**
	 * Creates a MTMCTS player that simulates random plays on separate threads
	 * and estimates the best position on the game grid with a given depth limit and a default verbosity.
	 * @param depthLimit The max depth limit for each simulated play
	 * @param threadMultiplier The number of threads run per turn (available positions multiplied by this)
	 */
	public MTMCTSPlayer(int depthLimit, int threadMultiplier)
	{
		this(depthLimit, threadMultiplier, DEFAULT_VERBOSITY);
	}

	/**
	 * Creates a MTMCTS player that simulates random plays on separate threads
	 * and estimates the best position on the game grid with a given depth limit and a given verbosity.
	 * @param verbose Flag to print more detailed information while performing simulations and calculating results
	 */
	public MTMCTSPlayer(boolean verbose)
	{
		this(DEFAULT_DEPTH_LIMIT, DEFAULT_THREAD_MULTIPLIER, verbose);
	}

	/**
	 * Creates a MTMCTS player that simulates random plays on separate threads
	 * and estimates the best position on the game grid with a given depth limit and a given verbosity.
	 * @param depthLimit The max depth limit for each simulated play
	 * @param threadMultiplier The number of threads run per turn (available positions multiplied by this)1
	 * @param verbose Flag to print more detailed information while performing simulations and calculating results
	 */
	public MTMCTSPlayer(int depthLimit, int threadMultiplier, boolean verbose)
	{
		this.depthLimit = Math.min(depthLimit, NUM_POS);
		this.threadMultiplier = threadMultiplier;
		this.verbose = verbose;

		if(verbose)
		{
			System.out.println("Multi-threaded Monte Carlo Tree Search player");
			System.out.println("Depth Limit: "+depthLimit);
			System.out.println("Thread Multiplier: "+threadMultiplier);
		}
	}

	/**
	 * @see PokerSquaresPlayer#setPointSystem(PokerSquaresPointSystem, long)
	 */
	@Override
	public void setPointSystem(PokerSquaresPointSystem system, long millis)
	{
		this.system = system;
	}

	/**
	 * @see PokerSquaresPlayer#init()
	 */
	@Override
	public void init()
	{
		currentRound = 0;
		totalPlayTime = 0;
		random = new Random();
		
		// Get an array of all cards in the deck
		cardsInDeck = new ArrayList<>(Arrays.asList(Card.getAllCards()));

		availablePositions = new ArrayList<int[]>();
		grid = new Card[SIZE][SIZE];

		for (int x=0; x<SIZE; x++)
		{
			for (int y=0; y<SIZE; y++)
			{
				grid[x][y] = null;
				availablePositions.add(new int[]{x,y});
			}
		}
	}

	/**
	 * @see PokerSquaresPlayer#getPlay(Card, long)
	 */
	public int[] getPlay(Card card, long millisRemaining)
	{
		long startTime = System.currentTimeMillis();
		cardsInDeck.remove(card);
		int[] playPosition = availablePositions.get(availablePositions.size()-1);

		// Calculated total time this round can use
		long expectedTimePerRound = (millisRemaining / (NUM_POS - currentRound));
		long timeBuffer = (NUM_POS * 6);
		long calculatedTimePerRound = expectedTimePerRound - timeBuffer;

		currentRound++;

		// Only perform threaded simulations when not on last round
		if(currentRound < NUM_POS)
		{
			ArrayList<SimulationThread> results = runThreads(card, calculatedTimePerRound);
			playPosition = getBestPosition(results);
		}

		// Play card at position on the grid
		grid[playPosition[0]][playPosition[1]] = card;
		availablePositions.remove(playPosition);

		// Only calculate and print values if verbose is true
		if(verbose)
		{
			System.out.println("Time alloted: "+calculatedTimePerRound+" ms");
			// Update and print time for this round
			long timeUsed = (System.currentTimeMillis()-startTime);
			System.out.println("Time used:    "+timeUsed+" ms");

			// Update play time and print total time on last round
			totalPlayTime += timeUsed;
			if(currentRound == NUM_POS)
			{
				System.out.println("Total time taken: "+totalPlayTime+" ms");
			}
		}

		return playPosition;
	}

	/**
     * Compile and run a list of threads based off of the number of remaining positions on the grid.
	 * @param card The card given to play for this turn
	 * @param timePerRound The amount of time given for the threads to run
	 * @return The resulting list of finished threads with their computed states about their given position
	 */
	private ArrayList<SimulationThread> runThreads(Card card, long timePerRound)
	{
		ExecutorService es = Executors.newCachedThreadPool();
		ArrayList<SimulationThread> threads = new ArrayList<SimulationThread>();

		List<Card> threadSafeCardsInDeck = Collections.synchronizedList(new ArrayList<Card>(cardsInDeck));
		List<int[]> threadSafeAvailablePositions = Collections.synchronizedList(new ArrayList<int[]>(availablePositions));
		
		// Start 2 thread simulations for each available position on the grid
		for(int i=0; i<availablePositions.size()*threadMultiplier; i++)
		{
			SimulationThread thread = new SimulationThread(grid, card, availablePositions.get(i/threadMultiplier), threadSafeCardsInDeck, threadSafeAvailablePositions, depthLimit);
			threads.add(thread);
			es.submit(thread);
		}
		es.shutdown();
		
		// Shutdown all threads after timeout
		try
		{
			es.awaitTermination(timePerRound, TimeUnit.MILLISECONDS);
		}
		catch(InterruptedException ex)
		{
			// Fail silently
		}
		finally
		{ 
			es.shutdownNow();
		}

		return threads;
	}

	/**
	 * Calculate the best result and play position based off of the results average score.
	 * @param results A list of threads that completed running that have computed information about each play position
	 * @return A valid x and y position for the game grid
	 */
	private int[] getBestPosition(ArrayList<SimulationThread> results)
	{
		// Result with the best heuristic score
		SimulationThread bestResult = null;
		// Total number of simulations this round, only used when verbose is true
		long totalSimulations = 0;

		for(SimulationThread result : results)
		{	
			// Only calculate statistics about each thread and current play if verbose
			if(verbose)
			{
				result.printMetrics();
				totalSimulations += result.totalSimulations;
			}

			if(bestResult == null || result.averageScore > bestResult.averageScore)
				bestResult = result;
		}
		// Only print statistics about each thread and current play if verbose
		if(verbose)
		{
			// Print information about this play
			System.out.println("Playing position: {"+bestResult.position[0]+","+bestResult.position[1]+"}");
			System.out.println("Best h(score)="+String.format("%.5f",bestResult.averageScore)+", sims="+bestResult.totalSimulations);
			System.out.println("Total simulations: "+totalSimulations);
		}

		return bestResult.position;
	}

	/**
	 * @see PokerSquaresPlayer#getName()
	 */
	@Override
	public String getName()
	{
		return "MTMCTS-Depth"+depthLimit+"-Threadx"+threadMultiplier;
	}

	public class SimulationThread implements Runnable
	{
		// This threads id
		public long threadId;

		// Grid to play cards, copy of global grid
		private Card[][] playGrid;
		// Card given to play this round
		private Card cardToPlay;
		// This threads position to play the given card
		public int[] position;

		// Copy of global cardsInDeck used for playing new cards at each depth
		private List<Card> remainingCards;
		// Copy of global availablePositions used for picking available positions at each depth
		private List<int[]> availablePositions;

		// Given depth limit that cannot be exceeded
		private int depthLimit;

		// Metrics
		// Total number times a full simulation is played (from depth 0 to depthLimit) 
		public long totalSimulations;
		// Total accumulated score of each simulation
		public long totalScore;
		// Average score after each simulation is run: score/simCount
		public double averageScore;

		/**
		 * Simulation threads are threads used for doing a Monte Carlo Tree Search simulation starting at a given position with a given card, then
		 * performing subsequent simulations with random positions and random cards.
		 * @param grid A copy of the game grid from the parent class, contains all cards at their already played positions
		 * @param card The card that must be plated for this round
		 * @param pos The position this thread is tasked with evaluating
		 * @param remainingCards A copy of the remaining cards from the parent class, contains all cards that have yet to be played
		 * @param availablePositions A copy of the available positions from the parent class, contains all positions that have yet to be played at
		 * @param depthLimit The maximum depth limit of the simulation
		 */
		public SimulationThread(Card[][] grid, Card card, int[] pos, List<Card> remainingCards, List<int[]> availablePositions, int depthLimit)
		{
			this.playGrid = Arrays.stream(grid).map(Card[]::clone).toArray(Card[][]::new);
			this.cardToPlay = card;
			this.position = pos;

			this.remainingCards = new ArrayList<Card>(remainingCards);
			this.availablePositions = new ArrayList<int[]>(availablePositions);
			this.availablePositions.remove(this.position);
			
			this.depthLimit = depthLimit;

			this.totalSimulations = 0;
			this.totalScore = 0;
			this.averageScore = 0;
		}

		/**
		 * Required method for classes that implement Runnable. Is called with a this class is submitted to 
		 * an ExecutorService.
		 */
		public void run()
		{
			// Get this threads id
			threadId = Thread.currentThread().getId();
			
			// Play initial card
			playGrid[position[0]][position[1]] = cardToPlay;
			
			// Run simulations for possible scores from playing at this threads given position
			try { simulatePlaysAtDepth(depthLimit); }
			catch(Exception ex) { ex.printStackTrace(); }

			// Calculate average score based on total accumulated score divided by number of simulations tried
			averageScore = (double)totalScore/(double)totalSimulations;
		}

		/**
		 * Performs a Monte Carlo Tree Search simulation at a given depth and calculate the total score as well as the 
		 * number of simulations performed. A simulation is a complete (or when interrupted, partial) run 
		 * from depth 0 to the max depth.
		 * @param depthLimit The depth limit that the simulation can will simulate to 
		 */
		public void simulatePlaysAtDepth(int depthLimit)
		{
			// Get the maximum depth limit based on how many positions there are left to be played
			// If there are fewer positions than the depth limit, update the limit to the number
			// of remaining positions
			int maxDepth = Math.min(depthLimit, availablePositions.size());

			// No need to run simulation, calculate grid score and exit
			if(maxDepth == 0)
			{
				totalScore += system.getScore(playGrid); 
				totalSimulations = 1;
				return;
			}

			// Create a set of synchronized thread-safe lists for tracking cards played 
			// at random positions at a given depth for a single simulation
			List<int[]> positionsPlayed = Collections.synchronizedList(new ArrayList<int[]>(maxDepth));
			List<Card> cardsPlayed = Collections.synchronizedList(new ArrayList<Card>(maxDepth));

			// Loop while the current thread has not been interrupted
			while(!Thread.interrupted())
			{
				// Update simulation tries
				totalSimulations++;

				// Loop through depth
				for(int d=0; d<maxDepth; d++)
				{
					// Check for thread interruption and exit here so not to waist time
					if(Thread.interrupted())
					{
						totalScore += system.getScore(playGrid);
						return;
					}

					// Select a random position and a random card
					int[] pos = selectRandom(availablePositions, positionsPlayed);
					Card card = selectRandom(remainingCards, cardsPlayed);
					// Play the card
					playGrid[pos[0]][pos[1]] = card;

				}
				
				// Update score
				totalScore += system.getScore(playGrid);

				// Undo previous plays
				for(int i=0; i<maxDepth; i++)
				{
					// Remove played cards from grid
					playGrid[positionsPlayed.get(i)[0]][positionsPlayed.get(i)[1]] = null;
					// Add played cards back to deck
					remainingCards.add(cardsPlayed.get(i));
					// Add position back to list of available positions
					availablePositions.add(positionsPlayed.get(i));
				}

				// Reset lists that track positions and cards for a simulation
				positionsPlayed.clear();
				cardsPlayed.clear();
			}
		}

		/**
		 * Randomly selects a value from a list and updates the passed in lists to maintin consistency
		 * through out each simulation.
		 * @param selectionList A list where the randomly selected value will be taken from
		 * @param seclectedList A list where the randomly selected value will be added to
		 * @return The randomly selected value
		 */
		private <T> T selectRandom(List<T> selectionList, List<T> seclectedList)
		{
			int index = random.nextInt(selectionList.size());
			T value = selectionList.get(index);

			seclectedList.add(value);
			selectionList.remove(index);
			
			return value;
		}
		
		/**
		 * Prints information about this thread.
		 * Note: Should only be called after the thread has completed it's run for the most
		 * accurate metrics.
		 */
		public void printMetrics()
		{
			String hscore = String.format("%.5f", (double)totalScore/totalSimulations);
			System.out.println("id="+threadId+", score="+totalScore+", sims="+totalSimulations+", h(score)="+hscore+", depth="+depthLimit+"\n");
		}
	}

	public static void main(String[] args)
	{
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getBritishPointSystem();
		System.out.println(system);
		new PokerSquares(new MTMCTSPlayer(14, 2, true), system).play();
	}
}