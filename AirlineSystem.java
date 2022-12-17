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
        Set<ArrayList<Route>>disResult = as.shortestDistanceItinerary(as.cityNames.get(0),as.cityNames.get(2));
        System.out.println(disResult);
        Set<ArrayList<Route>> priceResult = as.cheapestItinerary(as.cityNames.get(3),as.cityNames.get(8));
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
                path.add(new Route(cityNames.get(prev), cityNames.get(dest), G.distTo[dest], G.priceTo[dest]));
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
        HashSet<ArrayList<Route>> cheapestPaths = new HashSet<>();
        int s = cityNames.indexOf(source);
        int d = cityNames.indexOf(destination);
        int t = cityNames.indexOf(transit);

        // using the bfs method to find the shortest distance(weight)
        G.bfs(s, true, false);
        if(!G.marked[d]){
            return cheapestPaths;
        }else{
            cheapestPaths = flattenEdgeTo(s, t);
        }

        HashSet<ArrayList<Route>> cheapestPathsViaTransit = new HashSet<>();
        for(ArrayList<Route> path : cheapestPaths){
            for(Route route: path){
                if(route.source.equals(transit) || route.destination.equals(transit)){
                    cheapestPathsViaTransit.add(path);
                    break;
                }
            }
        }
        return cheapestPathsViaTransit;
    }

    public Set<Set<Route>> getMSTs(){
        // initializations
        Set<Set<Route>> MSTs = new HashSet<>();
        for(int i = 0; i < G.v; i++){
            if(!G.marked[i]){
                G.dijkstras(i, false);
            }
            MSTs.add(getMST(i));
        }

        return MSTs;
    }

    private Set<Route> getMST(int s) {
        Set<Route> MST = new HashSet<>();
        LinkedList<Integer> q = new LinkedList<>();
        q.add(s);

        while (!q.isEmpty()) {
            s = q.remove();
            for(int d = 0; d < G.v; d++){
                Set<Integer> edge = G.edgeTo[d];
                if(edge.contains(s)){
                    MST.add(new Route(cityNames.get(s), cityNames.get(d), G.distTo[d], G.priceTo[d]));
                    q.add(d);
                }
            }
        }
        return MST;
    }

    public Set<ArrayList<Route>> tripsWithin(String city, double budget) throws CityNotFoundException{

        Set<ArrayList<Route>> allPaths = new HashSet<>();
        if(!cityNames.contains(city))
            throw new CityNotFoundException(city);

        int s = cityNames.indexOf(city);
        G.bfs(s, true, false);
        for(int d = 0; d < G.v; d++){
            if(d != s && G.priceTo[d] <= budget){
                allPaths.addAll(flattenEdgeTo(s, d));
            }
        }

        return allPaths;
    }

    public Set<ArrayList<Route>> tripsWithin(double budget){
        Set<ArrayList<Route>> allPaths = new HashSet<>();

        for(int s = 0; s < G.v; s++){
            G.bfs(s, true, false);
            for(int d = 0; d < G.v; d++){
                if(d != s && G.priceTo[d] <= budget){
                    allPaths.addAll(flattenEdgeTo(s, d));
                }
            }
        }

        return allPaths;
    }

    private class Digraph {
        private final int v;
        private int e;
        private LinkedList<DirectedEdge>[] adj;
        private boolean[] marked;  // marked[v] = is there an s-v path
        private Set<Integer>[] edgeTo;      // edgeTo[v] = previous edges on shortest s-v path
        private int[] distTo;      // distTo[v] = weights of edges on shortest s-v path
        private double[] priceTo;      // priceTo[v] = prices of edges on shortest s-v path


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
            distTo = new int[this.v];
            edgeTo = new Set[this.v];
            priceTo = new double[this.v];

            Queue<Integer> q = new LinkedList<Integer>();
            for (int i = 0; i < v; i++){
                distTo[i] = INFINITY;
                marked[i] = false;
                priceTo[i] = Double.MAX_VALUE;
            }
            distTo[source] = 0;
            priceTo[source] = 0.0d;
            marked[source] = true;
            q.add(source);

            while (!q.isEmpty()) {
                int v = q.remove();
                for (DirectedEdge w : adj(v)) {
                    if (!marked[w.to()]) {
                        edgeTo[w.to()] = new HashSet<>();
                        edgeTo[w.to()].add(v);
                        if(useHopAsWeights){
                            distTo[w.to()] = distTo[v] + 1;
                        }else{
                            distTo[w.to()] = distTo[v] + w.weight;
                        }
                        priceTo[w.to()] = priceTo[v] + w.price;
                        marked[w.to()] = true;
                        q.add(w.to());
                    }else{
                        if(usePriceAsWeight){
                            if(priceTo[v] + w.price < priceTo[w.to()]){
                                edgeTo[w.to()] = new HashSet<>();
                                edgeTo[w.to()].add(v);
                                priceTo[w.to()] = priceTo[v] + w.price;
                                q.add(w.to());
                            }else if(priceTo[v] + w.price == priceTo[w.to()]){
                                edgeTo[w.to()].add(v);
                            }
                        }else if(useHopAsWeights){
                            if(distTo[v] + 1 < distTo[w.to()]){
                                edgeTo[w.to()] = new HashSet<>();
                                edgeTo[w.to()].add(v);
                                distTo[w.to()] = distTo[v] + 1;
                                q.add(w.to());
                            }else if(distTo[v] + 1 == distTo[w.to()]){
                                edgeTo[w.to()].add(v);
                            }
                        }else{
                            if(distTo[v] + w.weight < distTo[w.to()]){
                                edgeTo[w.to()] = new HashSet<>();
                                edgeTo[w.to()].add(v);
                                distTo[w.to()] = distTo[v] + w.weight;
                                q.add(w.to());
                            }else if(distTo[v] + w.weight == distTo[w.to()]){
                                edgeTo[w.to()].add(v);
                            }

                        }
                    }
                }
            }
        }

        public void dijkstras(int source, boolean usePriceAsWeight) {
            marked = new boolean[this.v];
            distTo = new int[this.v];
            edgeTo = new Set[this.v];
            priceTo = new double[this.v];


            for (int i = 0; i < v; i++){
                distTo[i] = INFINITY;
                marked[i] = false;
            }

            distTo[source] = 0;
            marked[source] = true;
            priceTo[source] = 0.0d;

            int nMarked = 1;

            int current = source;
            while (nMarked < this.v) {
                for (DirectedEdge w : adj(current)) {
                    if(!usePriceAsWeight){
                        if (distTo[current]+w.weight() < distTo[w.to()]) {
                            edgeTo[w.to()] = new HashSet<>();
                            edgeTo[w.to()].add(current);
                            distTo[w.to()] = distTo[current]+w.weight();
                        }
                    }else{
                        if (priceTo[current] + w.price() < priceTo[w.to()]) {
                            edgeTo[w.to()] = new HashSet<>();
                            edgeTo[w.to()].add(current);
                            priceTo[w.to()] = priceTo[current] + w.price();
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
                    if(!usePriceAsWeight && distTo[i] < min){
                        min = distTo[i];
                        current = i;
                    }
                    if(usePriceAsWeight && priceTo[i] < minDouble){
                        minDouble = priceTo[i];
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

