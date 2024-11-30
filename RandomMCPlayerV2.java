import java.util.ArrayList;
import java.util.Random;

/**
 * RandomMCPlayer - a simple Monte Carlo implementation of the player interface for PokerSquares.
 * For each possible play, continues play with random possible card draws and random card placements to a given depth limit 
 * (or game end).  Having sampled trajectories for all possible plays, the RandomMCPlayer then selects the
 * play yielding the best average scoring potential in such Monte Carlo simulation.
 * 
 * Disclaimer: This example code is not intended as a model of efficiency. (E.g., patterns from Knuth's Dancing Links
 * algorithm (DLX) can provide faster legal move list iteration/deletion/restoration.)  Rather, this example
 * code illustrates how a player could be constructed.  Note how time is simply managed so as to not run out the play clock.
 * 
 * Author: Todd W. Neller
 * Modifications by: Michael W. Fleming
 */
public class RandomMCPlayerV2 implements PokerSquaresPlayer {
	
	private final int SIZE = 5; // number of rows/columns in square grid
	private final int NUM_POS = SIZE * SIZE; // number of positions in square grid
	private final int NUM_CARDS = Card.NUM_CARDS; // number of cards in deck
	private Random random = new Random(); // pseudorandom number generator for Monte Carlo simulation 
	private int[] plays = new int[NUM_POS]; // positions of plays so far (index 0 through numPlays - 1) recorded as integers using row-major indices.
	// row-major indices: play (r, c) is recorded as a single integer r * SIZE + c (See http://en.wikipedia.org/wiki/Row-major_order)
	// From plays index [numPlays] onward, we maintain a list of yet unplayed positions.
	private int numPlays = 0; // number of Cards played into the grid so far
	private PokerSquaresPointSystem system; // point system
	private int depthLimit = 2; // default depth limit for Random Monte Carlo (MC) play
	private Card[][] grid = new Card[SIZE][SIZE]; // grid with Card objects or null (for empty positions)
	private Card[] simDeck = Card.getAllCards(); // a list of all Cards. As we learn the index of cards in the play deck,
	                                             // we swap each dealt card to its correct index.  Thus, from index numPlays 
												 // onward, we maintain a list of undealt cards for MC simulation.
	private int[][] legalPlayLists = new int[NUM_POS][NUM_POS]; // stores legal play lists indexed by numPlays (depth)
	// (This avoids constant allocation/deallocation of such lists during the selections of MC simulations.)

	/**
	 * Create a Random Monte Carlo player that simulates random play to depth 2.
	 */
	public RandomMCPlayerV2() {
	}
	
	/**
	 * Create a Random Monte Carlo player that simulates random play to a given depth limit.
	 * @param depthLimit depth limit for random simulated play
	 */
	public RandomMCPlayerV2(int depthLimit) {
		this.depthLimit = depthLimit;
	}
	
	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#init()
	 */
	@Override
	public void init() { 
		// clear grid
		for (int row = 0; row < SIZE; row++)
			for (int col = 0; col < SIZE; col++)
				grid[row][col] = null;
		// reset numPlays
		numPlays = 0;
		// (re)initialize list of play positions (row-major ordering)
		for (int i = 0; i < NUM_POS; i++)
			plays[i] = i;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#getPlay(Card, long)
	 */
	@Override
	public int[] getPlay(Card card, long millisRemaining) {
		/*
		 * With this algorithm, the player chooses the legal play that has the highest expected score outcome.
		 * This outcome is estimated as follows:
		 *   For each move, many simulated random plays to the set depthLimit are performed and the (sometimes
		 *     partially-filled) grid is scored.
		 *   For each play simulation, random undrawn cards are drawn in simulation and the player
		 *     picks a play position randomly.
		 *   After many such plays, the average score per simulated play is computed.  The play with the highest 
		 *     average score is chosen (breaking ties randomly).   
		 */
		
		// match simDeck to actual play event; in this way, all indices forward from the card contain a list of 
		//   undealt Cards in some permutation.
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;
		if (numPlays < 24) { // not the forced last play
			
			// compute average time per move evaluation
			int remainingPlays = NUM_POS - numPlays; // ignores triviality of last play to keep a conservative margin for game completion
			long millisPerPlay = millisRemaining / remainingPlays; // dividing time evenly with future getPlay() calls
			long millisPerMoveEval = millisPerPlay / remainingPlays; // dividing time evenly across moves now considered
			
			System.out.println("Round Time: "+millisPerPlay+" ms");
			
			// copy the play positions (row-major indices) that are empty
			System.arraycopy(plays, numPlays, legalPlayLists[numPlays], 0, remainingPlays);
			
			ArrayList<Integer> allPlays = new ArrayList<Integer>(remainingPlays);
			ArrayList<Double> allScores = new ArrayList<Double>(remainingPlays);
			long totalSims = 0;
			
			for (int i = 0; i < remainingPlays; i++) { // for each legal play position
				
				int play = legalPlayLists[numPlays][i];
				
				long startTime = System.currentTimeMillis();
				long endTime = startTime + millisPerMoveEval; // compute when MC simulations should end
				
				// Play original card at given position
				makePlay(card, play / SIZE, play % SIZE);
				
				int[] scoreAndSim = new int[2];
				simPlay(depthLimit, endTime, scoreAndSim);

				// update (if necessary) the maximum average score and the list of best plays
				double averageScore = (double) scoreAndSim[0] / scoreAndSim[1];
				totalSims += scoreAndSim[1];

				System.out.println("score="+scoreAndSim[0]+", sims="+scoreAndSim[1]+", h(score)="+String.format("%.5f", averageScore));
				allPlays.add(play);
				allScores.add(averageScore);

				undoPlay(); // undo the play under evaluation
			}

			Double maxAverageScore = Double.NEGATIVE_INFINITY;
			int bestPlay = 0;
			for(int i=0; i<allPlays.size(); i++)
			{
				if(allScores.get(i) > maxAverageScore)
				{
					maxAverageScore = allScores.get(i);
					bestPlay = allPlays.get(i);
				}
			}

			System.out.println("Best h(score)="+String.format("%.5f",maxAverageScore)+", totalSims="+totalSims);
			
			// update our list of plays, recording the chosen play in its sequential position; all onward from numPlays are empty positions
			int bestPlayIndex = numPlays;
			while (plays[bestPlayIndex] != bestPlay)
				bestPlayIndex++;
			plays[bestPlayIndex] = plays[numPlays];
			plays[numPlays] = bestPlay;
		}

		int[] playPos = {plays[numPlays] / SIZE, plays[numPlays] % SIZE}; // decode it into row and column
		makePlay(card, playPos[0], playPos[1]); // make the chosen play (not undoing this time)
		return playPos; // return the chosen play
	}

	/**
	 * From the chosen play, perform simulated Card draws and random placement (depthLimit) iterations forward 
	 * and return the resulting grid score.
	 * @param depthLimit - how many simulated random plays to perform
	 * @return resulting grid score after random MC simulation to given depthLimit
	 */
	private void simPlay(int depthLimit, long endTime, int[] scoreAndSim) {
		int simCount = 0;
		int scoreTotal = 0;
		while (System.currentTimeMillis() < endTime) { // perform as many MC simulations as possible through the allotted time
			// Perform a Monte Carlo simulation of random play to the depth limit or game end, whichever comes first.
			int depth = Math.min(depthLimit, NUM_POS - numPlays); // compute real depth limit, taking into account game end
			
			if (depth == 0) { // with zero depth limit, return current score
				scoreTotal += system.getScore(grid);
			}
			else {// up to the non-zero depth limit or to game end, iteratively make the given number of random plays 
				for (int d = 0; d < depth; d++) {

					// generate a random card draw
					Card card = simDeck[random.nextInt(NUM_CARDS - numPlays) + numPlays];
					
					// Get the legal plays
					int remainingPlays = NUM_POS - numPlays;
					System.arraycopy(plays, numPlays, legalPlayLists[numPlays], 0, remainingPlays);
					
					// choose a random play from the legal plays
					int play = legalPlayLists[numPlays][random.nextInt(remainingPlays)];
					
					// Play the card
					makePlay(card, play / SIZE, play % SIZE);
				}
				
				// Update Score
				scoreTotal += system.getScore(grid);
	
				// Undo MC plays.
				for (int d = 0; d < depth; d++) {
					undoPlay();
				}
			}	
			
			simCount++; // increment count of MC simulations
		}
		scoreAndSim[0] = scoreTotal;
		scoreAndSim[1] = simCount;
	}
	
	public void makePlay(Card card, int row, int col) {
		// match simDeck to event
		int cardIndex = numPlays;
		while (!card.equals(simDeck[cardIndex]))
			cardIndex++;
		simDeck[cardIndex] = simDeck[numPlays];
		simDeck[numPlays] = card;
		
		// update plays to reflect chosen play in sequence
		grid[row][col] = card;
		int play = row * SIZE + col;
		int j = 0;
		while (plays[j] != play)
			j++;
		plays[j] = plays[numPlays];
		plays[numPlays] = play;
		
		// increment the number of plays taken
		numPlays++;
	}

	public void undoPlay() { // undo the previous play
		numPlays--;
		int play = plays[numPlays];
		grid[play / SIZE][play % SIZE] = null;	
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#setPointSystem(PokerSquaresPointSystem, long)
	 */
	@Override
	public void setPointSystem(PokerSquaresPointSystem system, long millis) {
		this.system = system;
	}

	/* (non-Javadoc)
	 * @see PokerSquaresPlayer#getName()
	 */
	@Override
	public String getName() {
		return "RandomMCPlayerDepth" + depthLimit;
	}

	/**
	 * Demonstrate RandomMCPlay with Ameritish point system.
	 * @param args (not used)
	 */
	public static void main(String[] args) {
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getBritishPointSystem();
		System.out.println(system);
		new PokerSquares(new RandomMCPlayerV2(5), system).play(); // play a single game
	}

}
