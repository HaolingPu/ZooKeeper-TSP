//3E33912F8BAA7542FC4A1585D2DB6FE0312725B9
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <getopt.h>
#include <cmath>
#include <iomanip>
using namespace std;

enum class Mode {
    MST,
    FASTTSP,
    OPTTSP,
    Invalid
};

struct Cage {
    int x;
    int y;
    char type = 'b';
    bool visited = false;  //k
    double distance = numeric_limits<double>::infinity();   //d
    uint32_t parent = 0; //p
};


class Zookeeper {
    private:
    uint32_t countofb = 0;
    uint32_t countofs = 0;
    uint32_t countofw = 0;
    bool partb = false;
    double totalMSTWeight = 0; // for MST 
    double totalTSPWeight = 0;

    vector<Cage> cages;
    vector<uint32_t> tour;
    vector<vector<double>> distanceMatrix;


    // for permutation
    vector<uint32_t> bestpath;
    vector<uint32_t> currpath;
    double bestlength;
    double currlength = 0;





    public :
        Mode GetMode(int argc, char * argv[]) {
                int opt;
                Mode mode = Mode::Invalid;
                static struct option longOpts[] = {{ "mode", required_argument, nullptr, 'm' },
                                                { "help", no_argument, nullptr, 'h' },
                                                { nullptr, 0, nullptr, '\0' }};
                    while ((opt = getopt_long(argc, argv, "m:h", longOpts, nullptr)) != -1) {
                        switch (opt) {
                            case 'm':
                                   if (string(optarg) == "MST") {
                                    mode = Mode::MST;
                            
                                   }
                                    else if (string(optarg) == "FASTTSP") {
                                        mode = Mode::FASTTSP;
                                        partb = true;
                            
                                    }
                                    else if (string(optarg) == "OPTTSP") {
                                        mode = Mode::OPTTSP;
                                        
                                    }
                                    else {
                                        cerr << "Error: Invalid mode\n";
                                        exit(1);
                                    }
                                    break;
                            case 'h':
                                    cout << "I think you need some help!\n";  
                                    exit(0);
                            default:
                                   break;
                                 
                        }
                    }

                    if (mode == Mode::Invalid) {
                        cerr << "Error: No mode specified.\n";
                        exit(1);
                    }
            
            return mode;
        }


        void readin() {
            uint32_t V;
            cin >> V;

            cages.resize(V);
            for (uint32_t i = 0; i < V; ++i) {
                cin >> cages[i].x >> cages[i].y;
                if (cages[i].x > 0 || cages[i].y >0) {
                    cages[i].type = 's';   // safe animals
                    countofs++;
                }
                else if (cages[i].x < 0 && cages[i].y < 0){
                    cages[i].type = 'w';  // wild animals 
                    countofw++;
                }
                else {
                    cages[i].type = 'b';
                    countofb++;
                }
            }
     
        }

        double distance(const Cage& cage1, const Cage& cage2){
            double netX = cage1.x - cage2.x;
            double netY = cage1.y - cage2.y;
            return sqrt((netX)*(netX) + (netY)*(netY));
        }

        double distance_2(const Cage& cage1, const Cage& cage2){
            double netX = cage1.x - cage2.x;
            double netY = cage1.y - cage2.y;
            return (netX)*(netX) + (netY)*(netY);
        }


        bool valid_bridge (const Cage& cage1, const Cage& cage2){
            if (cage1.type != 'b' && cage2.type != 'b'){
                if (cage1.type != cage2.type) return false;
            }

            return true;
        }

    void primMST() {
            if (countofb == 0 && countofw != 0 && countofs != 0) {
                cerr << "Cannot construct MST\n" ;
                exit(1);
            }
            cages[0].distance = 0;
            
            for (uint32_t count = 0; count < cages.size(); ++count) {  
                double minWeight = numeric_limits<double>::infinity();
                uint32_t u = 0;

                for (uint32_t v = 0; v < cages.size(); ++v) {   // step1: loop through all and find the smallest false one
                    if (!cages[v].visited && cages[v].distance < minWeight) {
                        minWeight = cages[v].distance;
                        u = v;
                       
                    }
                }

                cages[u].visited = true; // step2: mark it as true

                totalMSTWeight += sqrt(minWeight);

                for (uint32_t v = 0; v < cages.size(); ++v) { // step3: update false neighbours (all nodes in same area)
                    if (!cages[v].visited && valid_bridge(cages[v], cages[u])) {
                        double weight = distance_2(cages[u], cages[v]);
                        if (weight < cages[v].distance) {
                            cages[v].parent = u;
                            cages[v].distance = weight;
                        }
                    }
                }
            }

                cout << totalMSTWeight << endl;
                for (uint32_t i = 0; i < cages.size(); i++){ 
                    // if (cages[i].parent != 0) {
                        if (i < cages[i].parent) cout << i << " " << cages[i].parent << endl;
                        else if (i > cages[i].parent) cout << cages[i].parent << " " << i << endl;
                    
                }
    }


double Insertionpath () {
    tour = {0,0};
    cages[0].visited = true;
    double netdistance;
        uint32_t nearestcage = 0;
        uint32_t insertplace = 0;
        for (uint32_t j = 1; j < cages.size(); j++){  // the node I am gon to insertm
        double minimalincrease = numeric_limits<double>::infinity();
            if (!cages[j].visited) {
                for (uint32_t i = 0; i < tour.size() - 1; i++) { // the every interval (i, i+1)  that i try to insert j.
                    netdistance = distance(cages[tour[i]], cages[j])
                    + distance(cages[j], cages[tour[i+1]])
                    - distance(cages[tour[i]], cages[tour[i+1]]);
                    if (netdistance < minimalincrease) {
                        minimalincrease = netdistance;
                        nearestcage = j;
                        insertplace = i + 1;
                    }
                }
        
            }
        tour.insert(tour.begin() + insertplace, nearestcage);
        cages[nearestcage].visited = true;
        }


    for (uint32_t i = 0; i < tour.size()-1 ;i++){
        totalTSPWeight += distance(cages[tour[i]], cages[tour[i+1]]);
    }
    tour.pop_back(); // 去掉最后一个0

    if (partb) {
        cout << totalTSPWeight << endl;
        for (uint32_t i = 0; i < tour.size(); i++){   // don't need to print the original 0.
            cout << tour[i] << " ";
        }
        cout << endl;
    }

    return totalTSPWeight;
}


void DistanceMatrix () {
    distanceMatrix.resize(cages.size(),vector<double>(cages.size(), 0.0));
    for (size_t i = 0; i < cages.size(); ++i) {
        for (size_t j = 0; j < cages.size(); ++j) {
            if (i != j) {
                distanceMatrix[i][j] = distance(cages[i], cages[j]);
            }
        }
    }
}

// promising 加一个if （k 《 5） promising 直接true
// k = path,size()- permlength
// use a functor for part A to have two version of MST
// 计算MST的时候，你的MST应该都是unvisited nodes， 可能会搞错的

double MST_new(uint32_t &permlength) {
        double totalweight = 0;
        size_t unvisits = currpath.size() - permlength;
        vector<bool> visited(unvisits, false);
        vector<double> distances(unvisits, numeric_limits<double>::infinity());
        distances[0] = 0;
        for (uint32_t i = 0; i < unvisits; i++){
            double minDistance = std::numeric_limits<double>::infinity();
            uint32_t vertex = 0;
            for (uint32_t j = 0; j < unvisits; j++) {
                if (!visited[j] && distances[j] < minDistance){
                    vertex = j;
                    minDistance = distances[j];
                }
            }

            visited[vertex] = true;
            totalweight += minDistance;

            for (uint32_t k =0; k < unvisits; k++) {
                if (!visited[k]) {
                    // double dis = distance_2(cages[currpath[permlength + vertex]], cages[currpath[permlength + k]]);
                    double dis = distanceMatrix[currpath[permlength + vertex]][currpath[permlength + k]];
                    if (dis < distances[k]){
                        distances[k] = dis;
                    }
                }
            }

        }     
        return totalweight;
    }


bool promising(uint32_t & permlength) {
   
    size_t k = currpath.size() - permlength;
    if (k < 5) return true;
    double mst = MST_new(permlength);
    double blueedge1 = numeric_limits<double>::infinity();
    double blueedge2 = numeric_limits<double>::infinity();
    for (uint32_t i = permlength; i < currpath.size() ; i++) {
        // double e1 = distance_2(cages[currpath[0]], cages[currpath[i]]);
        // double e2 = distance_2(cages[currpath[permlength-1]], cages[currpath[i]]);
        double e1 = distanceMatrix[currpath[0]][currpath[i]];
        double e2 = distanceMatrix[currpath[permlength-1]][currpath[i]];
        if (e1 < blueedge1) blueedge1 = e1;
        if (e2 < blueedge2) blueedge2 = e2;
    }
    // double estimate = sqrt(blueedge1) + sqrt(blueedge2) + mst + currlength;
    double estimate = blueedge1 + blueedge2 + mst + currlength;



//     for (size_t i = 0; i < currpath.size(); ++i){
//     cerr << setw(2) << currpath[i] << ' ';
// }
// cerr << setw(4) << permlength << setw(10) << currlength;
// cerr << setw(10) << blueedge1 << setw(10) << blueedge2;
// cerr << setw(10) << mst << setw(10) << estimate << "  " << (estimate < bestlength) << endl;
return estimate < bestlength;
}

void partC() {
    DistanceMatrix();
    bestlength = Insertionpath();
    bestpath = tour;
    currpath = tour;
    genPerms(1);
    cout << bestlength << '\n';
    for (uint32_t i = 0; i < bestpath.size(); i++) {
        cout << bestpath[i] << " ";
    }
    cout << endl;
}

void genPerms(uint32_t permlength) {
    if (permlength == currpath.size()) {
        double closeEdge = distance(cages[currpath.back()], cages[currpath[0]]);
        // double closeEdge = distanceMatrix[currpath.back()][currpath[0]];
        currlength += closeEdge;
        if (currlength < bestlength){
            bestpath = currpath;
            bestlength = currlength;
        }
        currlength -= closeEdge;
        return;
    }  // if ..complete path

    if (!promising(permlength)) {
        return;
    }  // if ..not promising

// curcost should stay at 0
  for (size_t i = permlength; i < currpath.size(); ++i) {
    std::swap(currpath[permlength], currpath[i]);
    currlength += distance(cages[currpath[permlength]], cages[currpath[permlength-1]]);
    // currlength += distanceMatrix[currpath[permlength]][currpath[permlength-1]];
    genPerms(permlength + 1);
    currlength -= distance(cages[currpath[permlength]], cages[currpath[permlength-1]]);
    // currlength -= distanceMatrix[currpath[permlength]][currpath[permlength-1]];
    std::swap(currpath[permlength], currpath[i]);
  }  // for ..unpermuted elements

}  // genPerms()




};

int main(int argc, char* argv[]) { 
    std::ios_base::sync_with_stdio(false);
    Zookeeper zoo;
    Mode mode = zoo.GetMode(argc, argv);
    cout << fixed << showpoint << setprecision(2) << boolalpha;
    
    switch (mode) {
        case Mode::MST:
            zoo.readin();
            zoo.primMST();
            break;
        case Mode::FASTTSP:
            zoo.readin();
            zoo.Insertionpath();
            break;
        case Mode::OPTTSP:
            zoo.readin();
            zoo.partC();
            break;
        case Mode::Invalid:
        std::cerr << "Invalid mode selected." << std::endl;
        exit(1);
    // cerr << fixed << showpoint << setprecision(2) << boolalpha;
    }

    return 0;
}