/**
 * Graph representation.
 *
 * @author Roger S. Zou
 * @version 2.0 07/29/2015
 */

import java.util.*;

public class JunctionGraph {

    private HashMap<Vertex, HashSet<Vertex>> adjList;  // implemented with adjacency list

    private void initialize() {
        this.adjList = new HashMap<Vertex, HashSet<Vertex>>();
    }

    /**
     * Constructor to create a graph with k edges specified by the vertices specified
     * in the source vertex and destination vertex matrix.
     * @param sourceVertices (k x 2) matrix of k source vertices. The i-th (1 x 2) matrix
     *                       indicates the source vertex coordinates for the i-th edge.
     *                       The corresponding destination vertex of the edge is the i-th
     *                       (1 x 2) matrix in destVertices.
     * @param destVertices (k x 2) matrix of k destination vertices, where the i-th entry
     *                     corresponds to the i-th entry of sourceVertices
     */

    public JunctionGraph(float[][] sourceVertices, float[][] destVertices) {

        initialize();

        // Validate inputs
        if (sourceVertices.length != destVertices.length) {
            throw new IllegalArgumentException("Input matrices MUST be the same size!");
        }

        // Add all vertices and edges to graph
        for (int i = 0; i < sourceVertices.length; i++) {
            Vertex v1 = new Vertex(sourceVertices[i]);
            Vertex v2 = new Vertex(destVertices[i]);
            this.addVertex(v1);
            this.addVertex(v2);
            this.addEdge(v1, v2);
        }

    }

    public ArrayList<int[][]> getAdjList(float[][] srcV, float[][] dstV) {
        int counter;

        // Add vertex key information to hashtable.
        Set<Vertex> vertexSet = this.Vset();
        HashMap<Vertex, Integer> indexV = new HashMap<Vertex, Integer>();
        counter = 1;
        for (Vertex v : vertexSet) {
            indexV.put(v, counter);
            counter++;
        }

        // Add the edge key information to hashtable.
        LinkedList<Edge> edgeSet = new LinkedList<Edge>();
        for (int i = 0; i < srcV.length; i++)
            edgeSet.add(new Edge(new Vertex(srcV[i]), new Vertex(dstV[i])));
        HashMap<Edge, LinkedList<Integer>> indexE = new HashMap<Edge, LinkedList<Integer>>();
        counter = 1;
        for (Edge e : edgeSet) {
            LinkedList<Integer> entry;
            if (indexE.containsKey(e))
                entry = indexE.get(e);
            else
                entry = new LinkedList<Integer>();
            entry.addFirst(counter);
            indexE.put(e, entry);
            counter++;
        }

        // Construct adjacency list.
        ArrayList<int[][]> adjList = new ArrayList<int[][]>();
        for (Vertex u : vertexSet) {
            Set<Vertex> neighbors = this.getEdges(u);
            ArrayList<int[]> allist = new ArrayList<int[]>();
            for (Vertex v : neighbors) {
                int vInd = indexV.get(v);
                int eInd;
                LinkedList<Integer> entry = indexE.get(new Edge(u,v));
                for (int j = 0; j < entry.size(); j++) {
                    eInd = entry.get(j);
                    allist.add(new int[]{vInd, eInd});
                }
            }
            int[][] adjlisti = new int[allist.size()][2];
            for (int i = 0; i < allist.size(); i++)
                adjlisti[i] = allist.get(i);
            adjList.add(adjlisti);
        }
        return adjList;
    }


    public ArrayList<int[][]> getFaceList(float[][] srcV, float[][] dstV) {
        int counter;

        // Add vertex key information to hashtable.
        Set<Vertex> vertexSet = this.Vset();
        HashMap<Vertex, Integer> indexV = new HashMap<Vertex, Integer>();
        counter = 1;
        for (Vertex v : vertexSet) {
            indexV.put(v, counter);
            counter++;
        }

        // Add the edge key information to hashtable.
        LinkedList<Edge> edgeSet = new LinkedList<Edge>();
        for (int i = 0; i < srcV.length; i++)
            edgeSet.add(new Edge(new Vertex(srcV[i]), new Vertex(dstV[i])));
        HashMap<Edge, LinkedList<Integer>> indexE = new HashMap<Edge, LinkedList<Integer>>();
        counter = 1;
        for (Edge e : edgeSet) {
            LinkedList<Integer> entry;
            if (indexE.containsKey(e))
                entry = indexE.get(e);
            else
                entry = new LinkedList<Integer>();
            entry.addFirst(counter);
            indexE.put(e, entry);
            counter++;
        }

        ArrayList<EN> enList = new ArrayList<EN>();
        for (Vertex root : vertexSet) {
            enList.add(new EN(root, this.getEdges(root)));
        }

        HashSet<Cell> allCells = new HashSet<Cell>();
        for (Edge e : edgeSet) {
            ArrayList<Vertex> vListN = new ArrayList<Vertex>();
            ArrayList<Edge> eListN = new ArrayList<Edge>();
            vListN.add(e.getLast());
            vListN.add(e.getFirst());
            eListN.add(e);
            getFaceListNextHelper(e.getLast(), e.getFirst(), enList, indexV, vListN, eListN);
            if (vListN.size() != eListN.size()) {
                throw new IllegalArgumentException("vListN != eListN");
            }
            allCells.add(new Cell(vListN, eListN));

            ArrayList<Vertex> vListB = new ArrayList<Vertex>();
            ArrayList<Edge> eListB = new ArrayList<Edge>();
            vListB.add(e.getLast());
            vListB.add(e.getFirst());
            eListB.add(e);
            getFaceListBackHelper(e.getLast(), e.getFirst(), enList, indexV, vListB, eListB);
            if (vListB.size() != eListB.size()) {
                throw new IllegalArgumentException("vListB != eListB");
            }
            allCells.add(new Cell(vListB, eListB));
        }

        ArrayList<int[][]> faceList = new ArrayList<int[][]>();
        for (Cell c : allCells) {
            ArrayList<Vertex> listV = c.getVList();
            ArrayList<Edge> listE = c.getEList();
            Iterator<Vertex> vIt = listV.iterator();
            Iterator<Edge> eIt = listE.iterator();
            ArrayList<int[]> faceListi = new ArrayList<int[]>();
            while (vIt.hasNext()) {
                int vKey = indexV.get(vIt.next());
                LinkedList<Integer> eKeys = indexE.get(eIt.next());
                for (Integer eKey : eKeys) {
                    faceListi.add(new int[]{vKey, eKey});
                }
            }
            int[][] faceListiarray = new int[faceListi.size()][2];
            for (int i = 0; i < faceListi.size(); i++) {
                faceListiarray[i] = faceListi.get(i);
            }
            faceList.add(faceListiarray);
        }
        return faceList;

    }

    private void getFaceListNextHelper(Vertex u, Vertex root, ArrayList<EN> enList,
                                   HashMap<Vertex, Integer> indexV,
                                   ArrayList<Vertex> vList, ArrayList<Edge> eList) {
        int rootIndex = indexV.get(root) - 1;
        EN en = enList.get(rootIndex);
        Vertex v = en.getNext(u);
        if (vList.get(0).equals(v)) {
            eList.add(new Edge(root, v));
        } else {
            vList.add(v);
            eList.add(new Edge(root, v));
            getFaceListNextHelper(root, v, enList, indexV,  vList, eList);
        }
    }

    private void getFaceListBackHelper(Vertex u, Vertex root, ArrayList<EN> enList,
                                       HashMap<Vertex, Integer> indexV,
                                       ArrayList<Vertex> vList, ArrayList<Edge> eList) {
        int rootIndex = indexV.get(root) - 1;
        EN en = enList.get(rootIndex);
        Vertex v = en.getBack(u);
        if (vList.get(0).equals(v)) {
            eList.add(new Edge(root, v));
        } else {
            vList.add(v);
            eList.add(new Edge(root, v));
            getFaceListBackHelper(root, v, enList, indexV,  vList, eList);
        }
    }

    /**
     * @return the # of edges
     */
    public int edgeCount() {
        return this.Eset().size();
    }

    /**
     * @return the # of vertices
     */
    public int vertexCount() {
        return adjList.size();
    }

    /**
     * @return the set of vertices as Vertex objects
     */
    private Set<Vertex> Vset() {
        return adjList.keySet();
    }

    /**
     * @return the (row, col) coordinates of vertices in a nx2 matrix.
     */
    public float[][] V() {
        float[][] vertexMatrix = new float[adjList.size()][2];
        int counter = 0;
        for (Vertex n : adjList.keySet()) {
            vertexMatrix[counter] = n.toArray();
            counter++;
        }
        return vertexMatrix;
    }

    /**
     * Finds all the edges in the undirected graph
     * @return the set of edges as Edge objects
     */
    private Set<Edge> Eset() {
        Set<Edge> edgeSet = new HashSet<Edge>();
        for (Vertex u : adjList.keySet()) {
            for (Vertex v: adjList.get(u)) {
                edgeSet.add(new Edge(u, v));
            }
        }
        return edgeSet;
    }

    /**
     * Checks if there exists an edge from u to v
     * @param u Vertex u
     * @param v Vertex v
     * @return true if there exists e(u,v). false otherwise.
     */
    public boolean hasEdge(Vertex u, Vertex v) {
        if (!adjList.containsKey(u) || !adjList.containsKey(v))
            return false;
        else
            return adjList.get(u).contains(v);
    }

    /**
     * Adds a vertex to the graph
     * @param v Vertex to be inserted into graph
     * @return true if there are no other vertices with the same row, col coordinate, false otherwise.
     */
    public boolean addVertex(Vertex v) {
        if (adjList.containsKey(v))
            return false;
        else {
            adjList.put(v, new HashSet<Vertex>());
            return true;
        }
    }


    /**
     * Update the position of a Vertex. NOTE that if the new position is identical to that of a vertex
     * already in the graph, then the two vertices are 'merged'.
     * @param vOld old Vertex
     * @param vNew new Vertex
     * @return true if vertex coordinates are successfully updated, false otherwise.
     */
    public boolean updateVertex(Vertex vOld, Vertex vNew) {
        if (!this.adjList.containsKey(vOld))
            return false;
        else {
            // Iterate over each key
            ArrayList<Vertex> keyList = new ArrayList<Vertex>(this.adjList.keySet());
            for (Vertex v : keyList) {
                HashSet<Vertex> edges = this.adjList.get(v);
                // Replace old (key, value) pair if key is the old vertex with new vertex.
                if (v.equals(vOld)) {
                    this.adjList.remove(vOld);
                    this.adjList.put(vNew, edges);
                } else if (edges.contains(vOld)) {
                    edges.remove(vOld);
                    edges.add(vNew);
                    this.adjList.remove(v);
                    this.adjList.put(v, edges);
                }
            }
            // Clean up graph so that vertices don't point to themselves.
            keyList = new ArrayList<Vertex>(this.adjList.keySet());
            for (Vertex v : keyList) {
                HashSet<Vertex> edges = this.adjList.get(v);
                if (edges.contains(v)) {
                    edges.remove(v);
                    this.adjList.remove(v);
                    this.adjList.put(v, edges);
                }
            }
        }

        return true;
    }

    /**
     * Update the position of a Vertex
     * @param rowOld row of old Vertex
     * @param colOld column of old Vertex
     * @param rowNew row of updated Vertex
     * @param colNew column of updated Vertex
     * @return true if vertex coordinates are successfully updated, false otherwise.
     */
    public boolean updateVertex(float rowOld, float colOld, float rowNew, float colNew) {
        return this.updateVertex(new Vertex(rowOld, colOld), new Vertex(rowNew, colNew));
    }

    /**
     * Returns all edges from a vertex
     * @param v Vertex v
     * @return Vertex u s.t there exists e(n,u).
     */
    public HashSet<Vertex> getEdges(Vertex v) {
        return adjList.containsKey(v) ? adjList.get(v) : null;
    }

    /**
     * Adds an edge between vertices.
     * @param u Vertex u
     * @param v Vertex v
     * @return true if a valid edge addition, false if an edge already exists.
     */
    public boolean addEdge(Vertex u, Vertex v) {
        if (hasEdge(u, v))
            return false;
        else {
            adjList.get(u).add(v);
            adjList.get(v).add(u);
            return true;
        }
    }

    /**
     * Remove an edge between vertices.
     * @param u Vertex u
     * @param v Vertex v
     * @return true if an existing edge between u and v is removed, false if there does not exist
     * an edge between u and v.
     */
    public boolean removeEdge(Vertex u, Vertex v) {
        if (!hasEdge(u, v)) return false;
        adjList.get(u).remove(v);
        adjList.get(v).remove(u);
        return true;
    }

    /**
     * Removes a vertex from graph
     * @param v Vertex to be removed from graph
     * @return true if vertex is successfully removed, false otherwise (probably didn't exist in graph)
     */
    public boolean removeVertex(Vertex v) {
        if (!adjList.keySet().contains(v))
            return false;
        else {
            adjList.remove(v);
            for (Vertex vKey : adjList.keySet()) {
                Set<Vertex> hashSet = adjList.get(vKey);
                if (hashSet.contains(v)) {
                    hashSet.remove(v);
                }
            }
            return true;
        }
    }

    /**
     * Encapsulation of a graph vertex.
     */
    class Vertex implements Comparable<Vertex> {

        private float row, col;
        private Vertex pointer;   // generic pointer useful for finding path in graph searching.

        public Vertex(float row, float col) {
            this.row = row;
            this.col = col;
            pointer = null;
        }

        public Vertex(float[] coords) {
            this.row = coords[0];
            this.col = coords[1];
            pointer = null;
        }

        public float getRow() { return row; }
        public float getCol() { return col; }

        public void set(float row, float col) {
            this.row = row;
            this.col = col;
        }

        public void clearPointer() { pointer = null; }
        public void setPointer(Vertex t) { pointer = t; }
        public Vertex getPointer() {return pointer; }

        /**
         * Computes Euclidean (L2) distance between 'this' vertex and another Vertex
         * @param v other vertex
         * @return euclidean distance
         */
        public double distL2(Vertex v) {
            return Math.sqrt( Math.pow(this.getRow()-v.getRow(),2) + Math.pow(this.getCol()-v.getCol(),2) );
        }

        /**
         * Convert Vertex to float array
         * @return float array of (row, col) coordinates
         */
        public float[] toArray() {
            return new float[] { this.row, this.col };
        }

        public int[] toIntArray() {
            return new int[] { Math.round(this.row), Math.round(this.col) };
        }

        @Override
        public int compareTo(Vertex b) {
            if (b == null) return 1;
            Vertex a = this;
            int aRow = (int) a.getRow();
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
            if (!(o instanceof Vertex)) return false;
            Vertex v = (Vertex) o;
            return (row == v.row && col == v.col);
        }

        @Override
        public int hashCode() {
            return (int) (row * 1187 + col * 2617);
        }

        @Override
        public String toString() {
            return "(" + this.getRow() + "," + this.getCol() + ")";
        }

    }

    /**
     * Encapsulation of a graph edge. Note that the ordering of the vertices DO NOT matter!
     * However, the first and last vertices in the edge are the most relevant in any case.
     */
    public class Edge {

        private LinkedList<Vertex> edge;

        public Edge(Vertex v1, Vertex v2) {
            edge = new LinkedList<Vertex>();
            if (v1.compareTo(v2) < 0) {
                edge.add(v2);
                edge.addFirst(v1);
            } else {
                edge.add(v1);
                edge.addFirst(v2);
            }
        }

        public Vertex getFirst() {
            return edge.getFirst();
        }
        public Vertex getLast() {
            return edge.getLast();
        }

        public LinkedList<Vertex> getEdge() {
            return edge;
        }
        public int getLength() {
            return edge.size();
        }

        @Override
        public boolean equals(Object o) {
            if (o == null) return false;
            if (o == this) return true;
            if (!(o instanceof Edge)) return false;
            Edge e = (Edge) o;

            // Edge must be same size and both sets have the same Vertex elements
            return new HashSet<Vertex>(this.getEdge()).equals(new HashSet<Vertex>(e.getEdge()));
        }

        @Override
        public int hashCode() {
            return this.getFirst().hashCode() + this.getLength() + this.getLast().hashCode();
        }

        @Override
        public String toString() {
            // Display only the first, last element, and length for brevity
            return "[" + this.getFirst().toString() + "->" + this.getLength() + "->" + this.getLast().toString() + "]";
        }

    }

    public class Cell {

        HashSet<Vertex> setV;
        HashSet<Edge> setE;
        ArrayList<Vertex> listV;
        ArrayList<Edge> listE;

        public Cell(ArrayList<Vertex> listV, ArrayList<Edge> listE) {
            this.setV = new HashSet<Vertex>(listV);
            this.setE = new HashSet<Edge>(listE);
            this.listV = listV;
            this.listE = listE;
        }

        public ArrayList<Vertex> getVList() { return this.listV; }
        public ArrayList<Edge> getEList() { return this.listE; }

        @Override
        public boolean equals(Object o) {
            if (o == null) return false;
            if (o == this) return true;
            if (!(o instanceof Cell)) return false;
            Cell c = (Cell) o;
            return (setV.equals(c.setV) && setE.equals(c.setE));
        }

        @Override
        public int hashCode() {
            int result = setV != null ? setV.hashCode() : 0;
            result = 31 * result + (setE != null ? setE.hashCode() : 0);
            return result;
        }
    }

    public class EN {

        Vertex root;
        ArrayList<Node> nodeList;
        ArrayList<Vertex> neighbors;

        public EN(Vertex root, Set<Vertex> neighbors) {
            this.root = root;
            this.neighbors = new ArrayList<Vertex>(neighbors);
            this.nodeList = new ArrayList<Node>();
            if (neighbors.size() == 1) {
                for (Vertex v : neighbors)
                    this.nodeList.add(new Node(root, v));
                Node n = this.nodeList.get(0);
                n.setNext(n);
                n.setBack(n);
            } else if (neighbors.size() > 1) {
                for (Vertex v : neighbors) {
                    this.nodeList.add(new Node(root, v));
                }
                Collections.sort(this.nodeList);
                for (int i = 0; i < this.nodeList.size()-1; i++) {
                    Node n1 = this.nodeList.get(i);
                    Node n2 = this.nodeList.get(i+1);
                    n1.setNext(n2);
                    n2.setBack(n1);
                }
                Node nbeg = this.nodeList.get(0);
                Node nend = this.nodeList.get(this.nodeList.size()-1);
                nbeg.setBack(nend);
                nend.setNext(nbeg);
            }
        }

        public Vertex getNext(Vertex v) {
            for (int i = 0; i < nodeList.size(); i++) {
                Node n = nodeList.get(i);
                if (v.equals(n.getVertex())) {
                    Node next = n.getNext();
                    if (next != null)
                        return next.getVertex();
                }
            }
            return null;
        }

        public Vertex getBack(Vertex v) {
            for (int i = 0; i < nodeList.size(); i++) {
                Node n = nodeList.get(i);
                if (v.equals(n.getVertex())) {
                    Node back = n.getBack();
                    if (back != null)
                        return back.getVertex();
                }
            }
            return null;
        }

        class Node implements Comparable<Node> {

            double theta;
            Vertex root, v;
            Node next, back;

            Node(Vertex root, Vertex v) {
                this.root = root;
                this.v = v;
                this.next = null;
                this.back = null;
                setTheta();
            }

            public double setTheta() {
                double ux = 1;
                double uy = 0;
                double vx = this.v.getCol() - this.root.getCol();
                double vy = this.v.getRow() - this.root.getRow();
                double norm = Math.sqrt(Math.pow(vx,2) + Math.pow(vy,2));
                vx = vx/norm;
                vy = vy/norm;
                this.theta = Math.acos(ux*vx + uy*vy);
                if (vy < 0) {
                    this.theta = 2*Math.PI - this.theta;
                }
                return this.theta;
            }

            public Vertex getVertex() {
                return v;
            }

            public Node getNext() {
                return next;
            }

            public Node getBack() {
                return back;
            }

            public void setNext(Node next) {
                this.next = next;
            }

            public void setBack(Node back) {
                this.back = back;
            }

            @Override
            public int compareTo(Node o) {
                double diff = this.theta - o.theta;
                if (diff < 0)
                    return -1;
                else if (diff > 0)
                    return 1;
                else
                    return 0;
            }

        }
    }

}
