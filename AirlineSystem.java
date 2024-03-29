import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.*;

class AirlineSystem implements AirlineInterface {
    private ArrayList<String> cityNames = null;
    private Digraph G = null;
    //private static Scanner scan = null;
    private static final int INFINITY = Integer.MAX_VALUE;
    // debugging purpose
    public static void main(String[] args) throws CityNotFoundException {
        AirlineSystem as = new AirlineSystem();
        String filename = "a4data1.txt";
        as.loadRoutes(filename);

//        System.out.println(as.cityNames);
//        HashSet<String> kk = new HashSet<>(as.retrieveCityNames());
//        System.out.println(kk);
//        int i=5;
//        System.out.println(as.cityNames.get(i));
//        HashSet<Route> route = new HashSet<>(as.retrieveDirectRoutesFrom(as.cityNames.get(i)));
//        System.out.println(route);
//        Set<ArrayList<String>> fewer = new HashSet<>(as.fewestStopsItinerary(as.cityNames.get(3),as.cityNames.get(8)));
//        System.out.println(fewer);
//        Set<ArrayList<Route>>disResult = as.shortestDistanceItinerary(as.cityNames.get(0),as.cityNames.get(2));
//        System.out.println(disResult);
//        Set<ArrayList<Route>> priceResult = as.cheapestItinerary(as.cityNames.get(3),as.cityNames.get(8));
//        System.out.println(priceResult);
//        Set<ArrayList<Route>> priceResult = as.cheapestItinerary(as.cityNames.get(3),as.cityNames.get(2), as.cityNames.get(8));
//        System.out.println(priceResult);
//        Set<ArrayList<Route>> priceResult = as.tripsWithin("Scranton",396);
//        System.out.println(priceResult);
//        Set<ArrayList<Route>> priceResult = as.tripsWithin(593.9733584102526);
//        System.out.println(priceResult);
        Set<Set<Route>> priceResult = as.getMSTs();
        System.out.println(priceResult);
    }

    public boolean loadRoutes(String fileName){
        try{
            Scanner fileScan = new Scanner(new FileInputStream(fileName));
            int v = Integer.parseInt(fileScan.nextLine());
            G = new Digraph(v);

            // using arraylist to store the cityNames
            cityNames = new ArrayList<>();
            for(int i=0; i<v; i++){
                cityNames.add(fileScan.nextLine());
            }

            while(fileScan.hasNext()){
                // store read values into from, to, weight, and price
                int from = fileScan.nextInt();
                int to = fileScan.nextInt();
                int weight = fileScan.nextInt();
                // price is double because it is a float type
                double price = fileScan.nextDouble();

                G.setDistTo(from - 1, to - 1, weight);
                G.setDistTo(to - 1, from - 1, weight);
                G.setPriceTo(from - 1, to - 1, price);
                G.setPriceTo(to - 1, from - 1, price);
                // use for debug
//                System.out.println(price);
//                System.out.println(weight);

                // added to the edge
                G.addEdge(new DirectedEdge(from - 1, to - 1, weight, price));
                G.addEdge(new DirectedEdge(to - 1, from - 1, weight, price));
            }
            return true;
        }catch(FileNotFoundException e){
            return false;
        }
    }

    public Set<String> retrieveCityNames(){
        return new HashSet<>(cityNames);
    }

    public Set<Route> retrieveDirectRoutesFrom(String city) throws CityNotFoundException{
        Set<Route> directRouts = new HashSet<>();
        // check if the request city is within the cityName list
        if(cityNames.contains(city)){
            // add to the route set -- directRouts
            for(DirectedEdge i: G.adj(cityNames.indexOf(city))){
                directRouts.add(new Route(city, cityNames.get(i.to()),i.weight, i.price));
                directRouts.add(new Route(cityNames.get(i.to()), city, i.weight, i.price));
            }
            // if not, throw exception
        }else
            throw new CityNotFoundException("Error");

        // return the result
        return directRouts;
    }

    public Set<ArrayList<String>> fewestStopsItinerary(String source, String destination) throws CityNotFoundException{
        // if source or destination are not in the list, throw exception
        if(!cityNames.contains(source))
            throw new CityNotFoundException(source);
        if (!cityNames.contains(destination))
            throw new CityNotFoundException(destination);
        // initializations
        HashSet<ArrayList<Route>> allPaths = new HashSet<>();
        HashSet<ArrayList<String>> result = new HashSet<>();
        int s = cityNames.indexOf(source);
        int d = cityNames.indexOf(destination);

        G.bfs(s, false, true);
        //System.out.println(s);

        // if no path to the request destination, return
        if(!G.marked[d]){
            return result;
        }else{
            // add the edges into cityPath list
            allPaths = flattenEdgeTo(s, d);
            for(ArrayList<Route> path: allPaths){
                ArrayList<String> pathString = new ArrayList<>();
                for(Route route: path){
                    pathString.add(route.destination);
                }
                pathString.add(0, source);
                result.add(pathString);
            }
        }
        return result;
    }

    public Set<ArrayList<Route>> shortestDistanceItinerary(String source, String destination) throws CityNotFoundException{
        if(!cityNames.contains(source))
            throw new CityNotFoundException(source);
        if (!cityNames.contains(destination))
            throw new CityNotFoundException(destination);

        // initializations
        HashSet<ArrayList<Route>> shortestDis = new HashSet<>();
        int s = cityNames.indexOf(source);
        int d = cityNames.indexOf(destination);

        // using the bfs method to find the shortest distance(weight)
        G.bfs(s, false, false);
        if(!G.marked[d]){
            return shortestDis;
        }else{
            shortestDis = flattenEdgeTo(s, d);
        }
        return shortestDis;
    }

    private HashSet<ArrayList<Route>> flattenEdgeTo(int source, int dest){
        if(source == dest){
            return new HashSet<>();
        }

        HashSet<ArrayList<Route>> thisRoutes = new HashSet<>();

        for(Integer prev: G.edgeTo[dest]){
            HashSet<ArrayList<Route>> prevRoutes = flattenEdgeTo(source, prev);
            if(prevRoutes.size() == 0){
                prevRoutes.add(new ArrayList<>());
            }
            for(ArrayList<Route> path: prevRoutes){
                path.add(new Route(cityNames.get(prev), cityNames.get(dest), G.distTo[prev][dest], G.priceTo[prev][dest]));
                thisRoutes.add(path);
            }
        }
        return thisRoutes;
    }
    public Set<ArrayList<Route>> cheapestItinerary(String source, String destination) throws CityNotFoundException{
        if(!cityNames.contains(source))
            throw new CityNotFoundException(source);
        if (!cityNames.contains(destination))
            throw new CityNotFoundException(destination);

        // initializations
        HashSet<ArrayList<Route>> cheapestPath = new HashSet<>();
        int s = cityNames.indexOf(source);
        int d = cityNames.indexOf(destination);

        // using the bfs method to find the shortest distance(weight)
        G.bfs(s, true, false);
        if(!G.marked[d]){
            return cheapestPath;
        }else{
            cheapestPath = flattenEdgeTo(s, d);
        }
        return cheapestPath;
    }

    public Set<ArrayList<Route>> cheapestItinerary(String source, String transit, String destination) throws CityNotFoundException{

        if(!cityNames.contains(source))
            throw new CityNotFoundException(source);
        if (!cityNames.contains(destination))
            throw new CityNotFoundException(destination);

        // initializations
        HashSet<ArrayList<Route>> allPaths = new HashSet<>();
        HashSet<ArrayList<Route>> transitPaths = new HashSet<>();
        int s = cityNames.indexOf(source);
        int d = cityNames.indexOf(destination);
        int t = cityNames.indexOf(transit);

        // using the bfs method to find the shortest distance(weight)
        G.bfs(s, true, false);
        if(!G.marked[t]){
            return allPaths;
        }else{
            transitPaths = flattenEdgeTo(s, t);
        }

        G.bfs(t, true, false);
        if(!G.marked[d]){
            return allPaths;
        }else{
            allPaths = flattenEdgeTo(t, d);
            for(ArrayList<Route> path : allPaths){

                for(ArrayList<Route> transitPath : transitPaths){
                    path.addAll(0, transitPath);
                }
            }
        }
        return allPaths;
    }

    public Set<Set<Route>> getMSTs(){
        // initializations to visit all nodes
        Set<Set<Route>> MSTs = new HashSet<>();
        boolean[] added = new boolean[G.v];
        for(int i = 0; i < G.v; i++){
            if(!added[i]){
                G.dijkstras(i, false);
                MSTs.add(getMST(i, added));
            }
        }

        return MSTs;
    }

    private Set<Route> getMST(int s, boolean[] added) {
        Set<Route> MST = new HashSet<>();
        LinkedList<Integer> q = new LinkedList<>();
        q.add(s);
        added[s] = true;

        while (!q.isEmpty()) {
            s = q.remove();
            for(int d = 0; d < G.v; d++){
                Set<Integer> edge = G.edgeTo[d];
                if(edge != null && edge.contains(s)){
                    MST.add(new Route(cityNames.get(s), cityNames.get(d), G.distTo[s][d], G.priceTo[s][d]));
                    q.add(d);
                    added[d] = true;
                }
            }
        }
        return MST;
    }

    public Set<ArrayList<Route>> tripsWithin(String city, double budget) throws CityNotFoundException{

        if(!cityNames.contains(city))
            throw new CityNotFoundException(city);

        int s = cityNames.indexOf(city);
        int mark[] = new int[G.v];
        Set<ArrayList<Route>> allPaths = tripsWithinRecursive(s, budget, mark);

        return allPaths;
    }

    private Set<ArrayList<Route>> tripsWithinRecursive(int s, double budget, int[] mark){
        Set<ArrayList<Route>> allPaths = new HashSet<>();
        if(budget <= 0){
            return allPaths;
        }


        mark[s] = 1;
        LinkedList<DirectedEdge> adjs = G.adj[s];
        for(DirectedEdge edge : adjs){
            int t = edge.to();
            if(mark[t] == 0 && budget - edge.price >= 0){
                for(ArrayList<Route> res : tripsWithinRecursive(t, budget - edge.price, mark)){
                    res.add(0, new Route(cityNames.get(s), cityNames.get(t), G.distTo[s][t], G.priceTo[s][t]));
                    allPaths.add(res);
                }

                // Default trip ending with t
                ArrayList<Route> path = new ArrayList<>();
                path.add(0, new Route(cityNames.get(s), cityNames.get(t), G.distTo[s][t], G.priceTo[s][t]));
                allPaths.add(path);
            }
        }

        mark[s] = 0;
        return allPaths;
    }

    public Set<ArrayList<Route>> tripsWithin(double budget){
        Set<ArrayList<Route>> allPaths = new HashSet<>();

        try {
            for(int s = 0; s < G.v; s++){
                allPaths.addAll(tripsWithin(cityNames.get(s), budget));
            }
        } catch (CityNotFoundException e) {
            e.printStackTrace();
        }

        return allPaths;
    }

    private class Digraph {
        private final int v;
        private int e;
        private LinkedList<DirectedEdge>[] adj;
        private boolean[] marked;  // marked[v] = is there an s-v path
        private Set<Integer>[] edgeTo;      // edgeTo[v] = previous edges on shortest s-v path
        private int[][] distTo;      // distTo[v] = weights of edges on shortest s-v path
        private double[][] priceTo;      // priceTo[v] = prices of edges on shortest s-v path


        /**
         * Create an empty digraph with v vertices.
         */
        public Digraph(int v) {
            if (v < 0) throw new RuntimeException("Number of vertices must be nonnegative");
            this.v = v;
            this.e = 0;
            @SuppressWarnings("unchecked")
            LinkedList<DirectedEdge>[] temp =
                    (LinkedList<DirectedEdge>[]) new LinkedList[v];
            adj = temp;
            for (int i = 0; i < v; i++)
                adj[i] = new LinkedList<DirectedEdge>();

            distTo = new int[this.v][this.v];
            priceTo = new double[this.v][this.v];
            for (int i = 0; i < v; i++){
                for(int j = 0; j < v; j++){
                    if(i == j){
                        distTo[i][j] = 0;
                        priceTo[i][j] = 0.0d;
                    }else{
                        distTo[i][j] = INFINITY;
                        priceTo[i][j] = Double.MAX_VALUE;
                    }
                }
            }
        }

        public void setDistTo(int s, int d, int dist){
            this.distTo[s][d] = dist;
        }

        public void setPriceTo(int s, int d, double price){
            priceTo[s][d] = price;
        }
        /**
         * Add the edge e to this digraph.
         */
        public void addEdge(DirectedEdge edge) {
            int from = edge.from();
            adj[from].add(edge);
            e++;
        }


        /**
         * Return the edges leaving vertex v as an Iterable.
         * To iterate over the edges leaving vertex v, use foreach notation:
         * <tt>for (DirectedEdge e : graph.adj(v))</tt>.
         */
        public Iterable<DirectedEdge> adj(int v) {
            return adj[v];
        }

        public void bfs(int source, boolean usePriceAsWeight, boolean useHopAsWeights) {
            marked = new boolean[this.v];
            edgeTo = new Set[this.v];

            Queue<Integer> q = new LinkedList<Integer>();
            marked[source] = true;
            q.add(source);

            while (!q.isEmpty()) {
                int v = q.remove();
                for (DirectedEdge w : adj(v)) {
                    if (!marked[w.to()]) {
                        edgeTo[w.to()] = new HashSet<>();
                        edgeTo[w.to()].add(v);
                        if(useHopAsWeights){
                            distTo[source][w.to()] = distTo[source][v] + 1;
                        }else{
                            distTo[source][w.to()] = distTo[source][v] + w.weight;
                        }
                        priceTo[source][w.to()] = priceTo[source][v] + w.price;
                        marked[w.to()] = true;
                        q.add(w.to());
                    }else{
                        if(usePriceAsWeight){
                            if(priceTo[source][v] + w.price < priceTo[source][w.to()]){
                                edgeTo[w.to()] = new HashSet<>();
                                edgeTo[w.to()].add(v);
                                priceTo[source][w.to()] = priceTo[source][v] + w.price;
                                q.add(w.to());
                            }else if(priceTo[source][v] + w.price == priceTo[source][w.to()]){
                                edgeTo[w.to()].add(v);
                            }
                        }else if(useHopAsWeights){
                            if(distTo[source][v] + 1 < distTo[source][w.to()]){
                                edgeTo[w.to()] = new HashSet<>();
                                edgeTo[w.to()].add(v);
                                distTo[source][w.to()] = distTo[source][v] + 1;
                                q.add(w.to());
                            }else if(distTo[source][v] + 1 == distTo[source][w.to()]){
                                edgeTo[w.to()].add(v);
                            }
                        }else{
                            if(distTo[source][v] + w.weight < distTo[source][w.to()]){
                                edgeTo[w.to()] = new HashSet<>();
                                edgeTo[w.to()].add(v);
                                distTo[source][w.to()] = distTo[source][v] + w.weight;
                                q.add(w.to());
                            }else if(distTo[source][v] + w.weight == distTo[source][w.to()]){
                                edgeTo[w.to()].add(v);
                            }

                        }
                    }
                }
            }
        }

        public void dijkstras(int source, boolean usePriceAsWeight) {
            marked = new boolean[this.v];
            edgeTo = new Set[this.v];
            int[] distToSource = new int[this.v];
            double[] priceToSource = new double[this.v];
            for (int i = 0; i < v; i++){
                if(i == source){
                    distToSource[i] = 0;
                    priceToSource[i] = 0.0d;
                }else{
                    distToSource[i] = INFINITY;
                    priceToSource[i] = Double.MAX_VALUE;
                }
            }

            marked[source] = true;
            int nMarked = 1;

            int current = source;
            while (nMarked < this.v) {
                for (DirectedEdge w : adj(current)) {
                    if(!usePriceAsWeight){
                        if (distToSource[current]+w.weight() < distToSource[w.to()]) {
                            edgeTo[w.to()] = new HashSet<>();
                            edgeTo[w.to()].add(current);
                            distToSource[w.to()] = distToSource[current]+w.weight();
                        }
                    }else{
                        if (priceToSource[current] + w.price() < priceToSource[w.to()]) {
                            edgeTo[w.to()] = new HashSet<>();
                            edgeTo[w.to()].add(current);
                            priceToSource[w.to()] = priceToSource[current] + w.price();
                        }
                    }
                }
                //Find the vertex with minimim path distance
                //This can be done more effiently using a priority queue!
                int min = INFINITY;
                double minDouble = Double.MAX_VALUE;
                current = -1;

                for(int i=0; i < this.v; i++){
                    if(marked[i])
                        continue;
                    if(!usePriceAsWeight && distToSource[i] < min){
                        min = distToSource[i];
                        current = i;
                    }
                    if(usePriceAsWeight && priceToSource[i] < minDouble){
                        minDouble = priceToSource[i];
                        current = i;
                    }
                }
                if(current >= 0){
                    marked[current] = true;
                    nMarked++;
                } else //graph is disconnected
                    break;
            }
        }
    }



    private static class DirectedEdge {
        private final int v;
        private final int w;
        private final int weight;
        private final double price;
        /**
         * Create a directed edge from v to w with given weight.
         */
        public DirectedEdge(int v, int w, int weight, double price) {
            this.v = v;
            this.w = w;
            this.weight = weight;
            this.price = price;
        }

        public int from(){
            return v;
        }

        public int to(){
            return w;
        }

        public int weight(){
            return weight;
        }

        public double price(){
            return price;
        }
    }
}


