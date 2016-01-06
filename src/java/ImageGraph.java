/**
 * Graph representation of an image with each pixel as a vertex, and each eight-neighborhood
 * connection as an edge.
 * Completely agnostic to vertex/edge encapsulations. Totally depends on matrices/arrays.
 * Supports standard graph operations.
 *
 * @author Roger S. Zou
 * @version 2.0 07/29/2015
 */

import java.util.*;

public class ImageGraph {

    private double[][] vertexCosts;
    private double[][][] edgeCosts;
    private boolean[][] validVertices;

    private void initialize(double[][] vertexCosts, double[][][] edgeCosts, boolean[][] validVertices) {
        // Initialize instance variables
        this.vertexCosts = vertexCosts;
        this.edgeCosts = edgeCosts;
        this.validVertices = validVertices;
    }

    /**
     * Constructor to create graph with each pixel as a vertex, with associated
     * vertex costs. Each vertex is connected to its 8 neighbors by an edge with
     * associated edge costs
     * @param vertexCosts matrix of vertex costs of size (m x n)
     * @param edgeCosts matrix of edge costs of size (m x n x 4). Only the bottom
     *                  left four vertices of 8 neighborhood is saved.
     * @param validVertices boolean matrix, where true indicates a vertex at a
     *                      particular index. All graph algorithms will only search
     *                      coordinates that are labeled true and treat them as
     *                      actual vertices.
     */
    public ImageGraph(double[][] vertexCosts, double[][][] edgeCosts, boolean[][] validVertices) {
        initialize(vertexCosts, edgeCosts, validVertices);
    }
    public ImageGraph(double[][] vertexCosts, double[][][] edgeCosts) {
        boolean[][] validVertices = new boolean[vertexCosts.length][vertexCosts[0].length];
        for (int i = 0; i < validVertices.length; i++) {
            for (int j = 0; j < validVertices[i].length; j++) {
                validVertices[i][j] = true;
            }
        }
        initialize(vertexCosts, edgeCosts, validVertices);
    }

    public int rowLength() {
        return this.vertexCosts.length;
    }
    public int colLength() {
        return this.vertexCosts[0].length;
    }

    /**
     * Check if a coordinate input is valid in the graph
     * @param row row of coordinate
     * @param col column of coordinate
     * @return true if the coordinates are within the bounds of 2D grid,
     * and a valid vertex as specified by validVertices
     */
    public boolean isValid(int row, int col) {
        return (0 <= row && row < this.rowLength() && 0 <= col &&
                col < this.colLength() && this.validVertices[row][col]);
    }

    /**
     * Shortest Path update function
     * @param vCost vertex cost
     * @param eCost edge cost
     * @return a function of the next vertex and edge costs
     */
    private double costFunctionUpdate(double vCost, double eCost) {
        return vCost + eCost;
    }

    /**
     * Compute the edge costs given source row and column, and list (row, col) of neighbors.
     * @param r1 source row index
     * @param c1 source column index
     * @param neighbors ArrayList of neighbor [row, col] indices
     * @return an array of edge costs whose index correspond to the index of neighbors list
     */
    private double[] getEdgeCosts(int r1, int c1, ArrayList<int[]> neighbors) {
        double[] eCosts = new double[neighbors.size()];
        int r2, c2;
        int counter = 0;
        for (int[] e : neighbors) {
            r2 = e[0];
            c2 = e[1];
            switch ((r2 - r1) * 10 + (c2 - c1)) {
                case 1 * 10 + 1:        eCosts[counter] = this.edgeCosts[r1][c1][3];  break;
                case 1 * 10 + 0:        eCosts[counter] = this.edgeCosts[r1][c1][2];  break;
                case 1 * 10 + (-1):     eCosts[counter] = this.edgeCosts[r1][c1][1];  break;
                case 0 * 10 + 1:        eCosts[counter] = this.edgeCosts[r1][c1][0];  break;
                case 0 * 10 + (-1):     eCosts[counter] = this.edgeCosts[r2][c2][0];  break;
                case (-1) * 10 + 1:     eCosts[counter] = this.edgeCosts[r2][c2][1];  break;
                case (-1) * 10 + 0:     eCosts[counter] = this.edgeCosts[r2][c2][2];  break;
                case (-1) * 10 + (-1):  eCosts[counter] = this.edgeCosts[r2][c2][3];  break;
            }
            counter++;
        }
        return eCosts;
    }

    /**
     * Implementation of Dijkstra's algorithm for single source shortest path with a min-queue
     * and early termination if all termination vertices are reached.
     * Time Complexity: O( |V| + |E|log|E| )
     * @param rowS source row index
     * @param colS source column index
     * @param rowT destination row indices
     * @param colT destination col indices
     * @return (m x n x 3) cost matrix C over entire graph. C(:,:,1) is the cost,
     * [ C(:,:,2), C(:,:,3) ] is the [ row, col ] of the previous vertex in path.
     */
    public double[][][] dijkstra(int rowS, int colS, int[] rowT, int[] colT) {

        // Initialize min queue Q and cost matrix C
        PriorityQueue<DijkstraNode> Q = new PriorityQueue<DijkstraNode>(8);
        double[][][] C = new double[this.rowLength()][this.colLength()][3];
        for (int r = 0; r < this.rowLength(); r++) {
            for (int c = 0; c < this.colLength(); c++) {
                if (r==rowS && c==colS) {
                    C[r][c][0] = this.vertexCosts[r][c];
                } else {
                    C[r][c][0] = Double.MAX_VALUE;
                }
                C[r][c][1] = -1;
                C[r][c][2] = -1;
            }
        }
        // Add source to queue
        Q.add(new DijkstraNode(this.vertexCosts[rowS][colS], rowS, colS));

        // Set termination conditions
        int terminationCount = rowT.length;

        // Loop while queue is not empty
        int r1, r2, c1, c2;
        double eCost, vCost, newCost;
        double[] eCosts;
        ArrayList<int[]> neighbors;
        while (!Q.isEmpty()) {
            DijkstraNode curVertex = Q.remove();
            r1 = curVertex.row;
            c1 = curVertex.col;
            neighbors = this.getNeighbors(r1, c1);
            eCosts = this.getEdgeCosts(r1, c1, neighbors);
            for (int e = 0; e < neighbors.size(); e++) {
                r2 = neighbors.get(e)[0];
                c2 = neighbors.get(e)[1];
                // Determine vertex and edge costs
                vCost = this.vertexCosts[r2][c2];
                eCost = eCosts[e];
                // calculate expected cost and update if necessary
                newCost = C[r1][c1][0] + this.costFunctionUpdate(vCost, eCost);
                if (newCost < C[r2][c2][0]) {
                    Q.remove(new DijkstraNode(C[r2][c2][0], r2, c2));
                    C[r2][c2][0] = newCost;
                    C[r2][c2][1] = r1;
                    C[r2][c2][2] = c1;
                    Q.add(new DijkstraNode(C[r2][c2][0], r2, c2));
                }
                // termination conditions
                for (int i = 0; i < rowT.length; i++) {
                    if (rowT[i]==r1 && colT[i]==c1) {
                        terminationCount--;
                        break;
                    }
                }
            }
        if (terminationCount < 0)
            break;
        }
        return C;
    }

    /**
     * Dijkstra's algorithm with single termination vertex
     * @param rowS source row index
     * @param colS source column index
     * @param rowT destination row index
     * @param colT destination col index
     * @return double[][][] dijkstra(rowS, colS, new int[]{ rowT }, new int[]{ colT });
     */
    public double[][][] dijkstra(int rowS, int colS, int rowT, int colT) {
        return this.dijkstra(rowS, colS, new int[]{ rowT }, new int[]{ colT });
    }

    /**
     * Class for use in Dijkstra's shortest path algorithm.
     * Used to determine node equality and node comparisons.
     * compareTo: two DijkstraNodes are comparable iff the cost are the same.
     * equals: two DijkstraNodes are equal iff the (row,col) coordinates are the same.
     * hashcode: if (row,col) coordinates are the same, then hash the same.
     */
    private class DijkstraNode implements Comparable<DijkstraNode> {

        private double cost;
        private int row, col;

        public DijkstraNode(double cost, int row, int col) {
            this.cost = cost;
            this.row = row;
            this.col = col;
        }


        @Override
        public int compareTo(DijkstraNode o) {
            if (this.cost < o.cost)
                return -1;
            else if (this.cost > o.cost)
                return 1;
            else
                return 0;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            DijkstraNode that = (DijkstraNode) o;

            return this.row == that.row && this.col == that.col;

        }

        @Override
        public int hashCode() {
            return 31 * this.row + 17 * this.col;
        }

        @Override
        public String toString() {
            return "(" + this.row + "," + this.col + ")";
        }
    }

    /**
     * Computes the path given the cost matrix with associated predecessor pointers generated
     * from shortest path algorithms
     * @param costMatrix a (m x n x 3) matrix, where mxn is the dimension of the image, and
     *                   (:,:,1) is the cost, (:,:,2) is the predecessor row, (:,:,3) is the
     *                   predecessor column.
     * @param rowT termination row index.
     * @param colT termination column index.
     * @return a (p x 2) matrix, where p is the length of the path from the source embedded
     * in costMatrix to the termination vertex, inclusive, and foreach i in p is a (1 x 2)
     * vector indicating the (row, col) of the specific point on the path.
     */
    public int[][] getPath(double[][][] costMatrix, int rowT, int colT) {
        LinkedList<int[]> path = new LinkedList<int[]>();

        // recursively find path from embedded source vertex to termination vertex
        getPathRecursiveHelper(costMatrix, rowT, colT, path);

        // convert from LinkedList<int[]> output to 2D matrix
        int[][] pathMatrix = new int[path.size()][2];
        int counter = 0;
        for (int[] aPath : path) {
            pathMatrix[counter] = aPath;
            counter++;
        }
        return pathMatrix;
    }

    /**
     * Finds the paths for all edges that minimizes cost function given edge inputs
     * @param startVertices (k x 2) matrix where the i-th (1 x 2) vector is
     *                      the start vertex of the i-th edge. There are k edges in total.
     *                      The end of the edge corresponds to the i-th (1 x 2) vector
     *                      in destVertices
     * @param destVertices (k x 2) matrix where the i-th (1 x 2) vector is the end
     *                     vertex of the i-th edge.
     * @return a size k ArrayList of (p x 2) matrices, for each path of size p, the (1 x 2)
     * vector corresponds to the (row, col) coordinate of the point on the path.
     */
    public ArrayList<int[][]> getEdgePaths(int[][] startVertices, int[][] destVertices) {

        int vCount = startVertices.length;
        ArrayList<int[][]> edgePaths = new ArrayList<int[][]>(vCount);

        int rowS, colS, rowT, colT;
        for (int i = 0; i < startVertices.length; i++) {
            rowS = startVertices[i][0];
            colS = startVertices[i][1];
            rowT = destVertices[i][0];
            colT = destVertices[i][1];
            edgePaths.add(this.getPath(this.dijkstra(rowS, colS, rowT, colT), rowT, colT));
        }
        return edgePaths;
    }

    /**
     * Recursive helper to getPath() that recursively computes the shortest path based on the cost matrix with
     * predecessor vertex information.
     * @param costMatrix a (m x n x 3) matrix, where mxn is the dimension of the image, and
     *                   (:,:,1) is the cost, (:,:,2) is the predecessor row, (:,:,3) is the
     *                   predecessor column.
     * @param curRow current row index.
     * @param curCol current column index.
     * @param path ArrayList storing the (row, col) coordinates of edges in the path.
     */
    public void getPathRecursiveHelper(double[][][] costMatrix, int curRow, int curCol, LinkedList<int[]> path) {
        path.addFirst(new int[]{ curRow, curCol });
        int nextRow = (int) costMatrix[curRow][curCol][1];
        int nextCol = (int) costMatrix[curRow][curCol][2];
        if (!(nextRow == -1 && nextCol == -1))
            this.getPathRecursiveHelper(costMatrix, nextRow, nextCol, path);
    }

    /**
     * Find the eight adjacent vertices to input vertex.
     * @param rowS source row index
     * @param colS source column index
     * @return an ArrayList of (row, col) points adjacent to input vertex.
     * Note that possibly less than 8 points found for corner cases.
     */
    public ArrayList<int[]> getNeighbors(int rowS, int colS) {
        ArrayList<int[]> neighbors = new ArrayList<int[]>();
        for (int r = rowS - 1; r <= rowS + 1; r++) {
            for (int c = colS - 1; c <= colS + 1; c++) {
                if ((r != rowS || c != colS) && isValid(r, c))
                    neighbors.add(new int[]{r, c});
            }
        }
        return neighbors;
    }

}
