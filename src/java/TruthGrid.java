/**
 * Generates edge information from ground truth data.
 *
 * @author Roger Zou
 * @version 2.0 07/29/2015
 */

import java.util.*;

public class TruthGrid {

    //constructor
    int[][] grid;

    public TruthGrid(int[][] imageMatrix) {
        this.grid = imageMatrix;
    }

    public int getNumRows() {
        return grid.length;
    }

    public int getNumCols() {
        return grid[0].length;
    }

    // returns true if the coordinates are within the bounds of 2D grid
    public boolean isValid(int row, int col) {
        return (0 <= row && row < getNumRows() && 0 <= col && col < getNumCols());
    }

    // returns true if coordinates are within bounds and grid is not empty at that location
    public boolean contains(int row, int col) {
        return (isValid(row, col) && grid[row][col] == 1);
    }

    //gets neighbors that are either one row OR one column away (no diagonals)
    public ArrayList<Point> getCrossNeighbors(int row, int col) {
        ArrayList<Point> neighbors = new ArrayList<Point>();
        if (contains(row - 1, col)) neighbors.add(new Point(row - 1, col));
        if (contains(row + 1, col)) neighbors.add(new Point(row + 1, col));
        if (contains(row, col - 1)) neighbors.add(new Point(row, col - 1));
        if (contains(row, col + 1)) neighbors.add(new Point(row, col + 1));
        return neighbors;
    }

    //gets neighbors that are only diagonal from center coordinate
    public ArrayList<Point> getDiagonalNeighbors(int row, int col) {
        ArrayList<Point> neighbors = new ArrayList<Point>();
        if (contains(row - 1, col - 1)) neighbors.add(new Point(row - 1, col - 1));
        if (contains(row - 1, col + 1)) neighbors.add(new Point(row - 1, col + 1));
        if (contains(row + 1, col - 1)) neighbors.add(new Point(row + 1, col - 1));
        if (contains(row + 1, col + 1)) neighbors.add(new Point(row + 1, col + 1));
        return neighbors;
    }

    // helper to check if a certain element in the grid is a node
    private boolean isJunction(int row, int col) {
        return findBranch(row, col).size() >= 3;
    }

    // finds all the vertices adjacent to a node that can extend as an Edge
    public ArrayList<Point> findBranch(int row, int col) {
        ArrayList<Point> branches = new ArrayList<Point>();
        if (contains(row, col)) {
            // cross neighbors must be branches
            ArrayList<Point> cn = getCrossNeighbors(row, col);
            for (Point v : cn)
                branches.add(v);
            // diagonal neighbors may be branches
            ArrayList<Point> dn = getDiagonalNeighbors(row, col);
            for (Point v : dn)
                if (!contains(v.getRow(), col) && !contains(row, v.getCol()))
                    branches.add(v);
        }
        return branches;
    }

    // finds the nodes (>=3 branches from Vertex) on the grid
    public ArrayList<Point> findJunctions() {
        ArrayList<Point> nodes = new ArrayList<Point>();
        for (int i = 0; i < getNumRows(); i++)
            for (int j = 0; j < getNumCols(); j++)
                if (contains(i, j) && isJunction(i, j))
                    nodes.add(new Point(i, j));

        return nodes;
    }

    // finds all junctions connected to this junction
    public ArrayList<Point> findEdges(Point node) {

        ArrayList<Point> edges = new ArrayList<Point>();
        ArrayList<Point> branches = findBranch(node.getRow(), node.getCol());    //the branches from center

        for (Point b : branches) {    //for each branch from node, make an edge segment

                LinkedList<Point> vertices = new LinkedList<Point>();
                HashSet<Point> visited = new HashSet<Point>();
                vertices.add(node);
                visited.add(node);

                if (findEdgeRecursive(b, vertices, visited))
                    edges.add(vertices.getLast());

        }
        return edges;
    }

    // forms an edge recursively, by adding to 'connector'. All coordinates traveled placed into visited.
    // returns true if a Connector is formed. returns false if a Fragment is formed
    private boolean findEdgeRecursive(Point v, List<Point> connector, Set<Point> visited) {

        // basis failure
        if (!contains(v.getRow(), v.getCol())) return false;    // runs out of bounds or no Vertex object exists

        // adds vertex to connector
        connector.add(v);
        visited.add(v);

        // gets the branches from v
        ArrayList<Point> branches = findBranch(v.getRow(), v.getCol());
        // basis failure: dead end, so this is a Fragment
        if (branches.size() <= 1) return false;

        // basis success: hits another node
        if (!v.equals(connector.get(0)) && isJunction(v.getRow(), v.getCol())) return true;

        for (Point b : branches) {
            if (!visited.contains(b)) {
                return findEdgeRecursive(b, connector, visited);
            }
        }
        return false;
    }

    /**
     * Compute all edges from image
     * @return a (2 x k x 2) matrices in an ArrayList, where the first matrix is (1,:,:), the second is
     * (2,:,:). The i-th component of the first matrix is a (1 x 2) vector of source vertex coordinates.
     * The second matrix contains the vertex coordinates for the destination vertices.
     */
    public int[][][] getAllEdges() {
        ArrayList<Point[]> allEdges = new ArrayList<Point[]>();

        ArrayList<Point> junctions = findJunctions();
        for (Point v : junctions) {
            ArrayList<Point> edges = findEdges(v);
            for (Point e : edges) {
                allEdges.add( new Point[]{v, e} );
            }
        }

        // convert from ArrayList of vertices to matrix
        int[][] sourceEdges = new int[allEdges.size()][2];
        int[][] destEdges = new int[allEdges.size()][2];
        for (int i = 0; i < allEdges.size(); i++) {
            Point[] edge = allEdges.get(i);
            sourceEdges[i] = new int[]{ edge[0].getRow(), edge[0].getCol() };
            destEdges[i] = new int[]{ edge[1].getRow(), edge[1].getCol() };
        }
        int[][][] allEdgesMatrix = new int[2][allEdges.size()][2];
        allEdgesMatrix[0] = sourceEdges;
        allEdgesMatrix[1] = destEdges;
        return allEdgesMatrix;
    }

    public class Point implements Comparable<Point> {

        //private variables
        private int row;
        private int col;

        public Point(int row, int col) {
            this.row = row;
            this.col = col;
        }

        //getters
        public int getRow() {
            return row;
        }

        public int getCol() {
            return col;
        }

        //Overridden methods
        @Override
        public int compareTo(Point b) {
            if (b == null) return 1;
            Point a = this;
            if (a.getRow() < b.getRow()) return -1;
            if (a.getRow() > b.getRow()) return 1;
            else {
                if (a.getCol() < b.getCol()) return -1;
                if (a.getCol() > b.getCol()) return 1;

                else return 0;
            }
        }

        @Override
        public boolean equals(Object o) {
            if (o == null) return false;
            if (o == this) return true;
            if (!(o instanceof Point)) return false;
            Point v = (Point) o;
            return (row == v.row && col == v.col);
        }

        @Override
        public int hashCode() {
            return row * 1187 + col * 2617;
        }

        @Override
        public String toString() {
            return "(" + this.getRow() + "," + this.getCol() + ")";
        }

        public void toPrint() {
            System.out.print("(" + this.getRow() + "," + this.getCol() + ")");
        }

    }
}
