#undef seed
#undef read
#undef write
#undef Move
#undef do_open
#undef do_close
#undef Null

#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <ctime>
#include <ratio>
#include <chrono>
#include <string>
#include <unistd.h>
#include <cstdio>
#include <stdio.h>
#include <getopt.h>
#include <experimental/optional>
#include <typeinfo>
#include "math.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <queue>
#include <cmath>
#include <spatialindex/capi/sidx_api.h>
#include <spatialindex/capi/sidx_impl.h>
#include <spatialindex/capi/sidx_config.h>
#include <ext/pb_ds/priority_queue.hpp>

#include <google/dense_hash_map>

#define LIMIT 0

using namespace std;
using namespace std::experimental;
using namespace std::chrono;
using namespace google;
using namespace SpatialIndex;

struct TreeSystem {
    double x;
    double y;
    double z;
    unsigned long int id64;
    bool isNeutron;
    bool isScoopable;
};

struct System {
    double x;
    double y;
    double z;
    unsigned long int id64;
    bool isNeutron;
    bool isScoopable;
    double distanceToDestination;
};

struct ShipDetails {
    double fuelMultiplier;
    double fuelPower;
    double baseMass;
    double optimalMass;
    double tankSize;
    double maxFuelPerJump;
    double rangeBoost;
    double scoopSpeed;
    bool injectionBoost;
};

struct Ship {
    double remainingFuel;
    bool isSupercharged;
    ShipDetails* base;
};

struct Node {
    System system;
    Node* parent;
    double f;
    double g;
    double h;
    Ship ship;
    double distance;
};

typedef dense_hash_map<unsigned long int,Node *> VisitedNode;
typedef dense_hash_map<Node *,bool> DeletedNode;

static inline double calculate_my_distance_squared(double &x1,double &y1,double &z1,double &x2,double &y2,double &z2) {
    return (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
}

static inline double calculate_my_distance(double &x1,double &y1,double &z1,double &x2,double &y2,double &z2) {
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
}

static inline double calculate_my_distance(const System& source,const System& destination) {
    double dx = destination.x-source.x;
    double dy = destination.y-source.y;
    double dz = destination.z-source.z;
    return sqrt(dx*dx+dy*dy+dz*dz);
}

static inline double calculate_my_distance_squared(const System& source,const System& destination) {
    double dx = destination.x-source.x;
    double dy = destination.y-source.y;
    double dz = destination.z-source.z;
    return dx*dx+dy*dy+dz*dz;
}

class SystemFetcher {
    protected:
        ISpatialIndex* tree;
        IStorageManager* diskfile;
        StorageManager::IBuffer* file;
    public:
        SystemFetcher() {
            char* pszVersion = SIDX_Version();
            cerr << "libspatialindex version: " << pszVersion << endl;
            free(pszVersion);

            string filename = "../GalaxySpatial/spatial_index";
            this->diskfile = StorageManager::loadDiskStorageManager(filename);
            this->file = StorageManager::createNewRandomEvictionsBuffer(*diskfile, 1024, false);

            this->tree = RTree::loadRTree(*file,1);

        }
        vector<System> lookup(System& source,System& destination,double range) {
            double plow[] = {source.x-range,source.y-range,source.z-range};
            double phigh[] = {source.x+range,source.y+range,source.z+range};
            // cerr << source.id64 << "(" << source.x << "," << source.y << "," << source.z << ")" << endl;
            Region r = Region(plow, phigh, 3);
            ObjVisitor* visitor = new ObjVisitor;
            tree->intersectsWithQuery(r, *visitor);
            std::vector<SpatialIndex::IData*>& results = visitor->GetResults();
            std::vector<System> capturedResults;
            double rangeSquared = range * range;
            for (auto &result : results) {
                uint32_t len;
                uint8_t* data;
                result->getData(len,&data);
                if (len == sizeof(TreeSystem)) {
                    TreeSystem *treeSystem = reinterpret_cast<TreeSystem*>(data);
                    System system;
                    system.id64 = treeSystem->id64;
                    system.x = treeSystem->x;
                    system.y = treeSystem->y;
                    system.z = treeSystem->z;
                    system.isNeutron = treeSystem->isNeutron;
                    system.isScoopable = treeSystem->isScoopable;
                    system.distanceToDestination = calculate_my_distance(system,destination);
                    double distanceSquared = calculate_my_distance_squared(source,system);
                    if (distanceSquared < rangeSquared) {
                        // cerr << system.id64 << "(" << system.x << "," << system.y << "," << system.z << ")" << endl;
                        capturedResults.push_back(system);
                    } else {
                        if (system.isNeutron) {
                            // cerr << "Skipping neutron at distance " << system.distanceToDestination << endl;
                        }
                        // cerr << "Distance was " << sqrt(distanceSquared) << " and should be less than " << sqrt(rangeSquared) << endl;
                    }
                } else {
                    cerr << "Length was " << len << " and should be " << sizeof(System) << endl;
                }
                delete[] data;
            }
            delete visitor;

            sort(capturedResults.begin(), capturedResults.end(), 
                [](const System& a, const System& b) -> bool { 
                    return a.distanceToDestination < b.distanceToDestination; 
                }
            );

            return capturedResults;
        }
        ~SystemFetcher() {
            delete this->tree;
            delete this->file;
            delete this->diskfile;
        }
};

double boostedFuelMultiplier(Ship ship) {
    double maxFuel = min(ship.remainingFuel,ship.base->maxFuelPerJump);
    double fuelMultiplier = ship.base->fuelMultiplier;
    double range = (ship.base->optimalMass / (ship.base->baseMass + ship.remainingFuel)) * (pow(maxFuel/fuelMultiplier, (1 / ship.base->fuelPower)));
    if (ship.isSupercharged) {
        range = range * 4;
    } else if (ship.base->injectionBoost) {
        range = range * 2;
    }
    if (ship.base->rangeBoost > 0) {
        fuelMultiplier = fuelMultiplier * pow(range / (range + ship.base->rangeBoost), ship.base->fuelPower);
    }
    return fuelMultiplier;
}

double calculateMaxRange(Ship ship) {
    double maxFuel = min(ship.remainingFuel,ship.base->maxFuelPerJump);
    double fuelMultiplier = boostedFuelMultiplier(ship);
    double range = (ship.base->optimalMass / (ship.base->baseMass + ship.remainingFuel)) * (pow(maxFuel/fuelMultiplier, (1 / ship.base->fuelPower)));
    if (ship.isSupercharged) {
        return range * 4;
    } else if (ship.base->injectionBoost) {
        range = range * 2;
    }
    return range;
}

double calculateMaxLadenRange(Ship ship) {
    double fuelMultiplier = boostedFuelMultiplier(ship);
    double range = (ship.base->optimalMass / (ship.base->baseMass + ship.base->tankSize)) * (pow(ship.base->maxFuelPerJump/fuelMultiplier, (1 / ship.base->fuelPower)));
    if (ship.isSupercharged) {
        return range*4;
    } else if (ship.base->injectionBoost) {
        range = range * 2;
    }
    return range;
}

static inline double calculate_cost_to_destination(const System& current,const Ship& ship) {
    double range = calculateMaxRange(ship);
    double ladenRange = range;
    if (ship.remainingFuel > ship.base->maxFuelPerJump) {
        ladenRange = calculateMaxLadenRange(ship);
    }

    double jumps = 1 + (current.distanceToDestination - range)/ladenRange;

    double refuel = 0;
    if (!current.isScoopable) {
        refuel += ((ship.base->tankSize - ship.remainingFuel) / ship.base->tankSize);
    }
    double cost = (jumps + refuel);
    /*
    if (ship.isSupercharged) {
        cerr << "Neutron found at cost " << cost << "(Jumps: " << jumps << " Refuel: " << refuel << ")" << endl;
        cerr << ship.base->tankSize << "," << ship.remainingFuel << "," << ship.base->maxFuelPerJump << "," << jumps << "," << ((ship.base->tankSize - ship.remainingFuel) / ship.base->maxFuelPerJump) << endl;
    }
    if (current.isScoopable) {
        cerr << "Scoopable with estimated cost " << cost << " at distance " << current.distanceToDestination << endl;
    } else {
        cerr << "Non Scoopable with estimated cost " << cost << " at distance " << current.distanceToDestination << endl;
    }
    */
    return cost;
}

inline double calculateFuelUsed(const Ship& ship,double distance) {
    double fuelMultiplier = boostedFuelMultiplier(ship);
    if (ship.isSupercharged) {
        distance = distance / 4;
    } else if (ship.base->injectionBoost) {
        distance = distance / 2;
    }
    return fuelMultiplier * pow(distance * ((ship.base->baseMass + ship.remainingFuel) / ship.base->optimalMass),ship.base->fuelPower);
}

inline double calculateRefuelingTime(const Ship& ship,double distance) {
    double fuelUsed = calculateFuelUsed(ship,distance);
    return (fuelUsed/ship.base->tankSize);
}

inline double calculateRefuelingTime(const Ship& ship) {
    return ((ship.base->tankSize - ship.remainingFuel)/ship.base->tankSize);
}

inline double calculateRefuelingTime(double tankSize,double fuelUsed) {
    return (fuelUsed/tankSize);
}

double averageJumpSize(Node* node) {
    unsigned long int count = 1;

    Node* current = node;
    double distance = node->system.distanceToDestination;
    while (current->parent != NULL) {
        count++;
        current = current->parent;
    }

    double result = (current->system.distanceToDestination - distance)/count;
    if (result < 0) {
        result = 0;
    }

    return result;
}

unsigned long int routeSize(Node* node) {
    unsigned long int count = 0;

    Node* current = node;
    while (current->parent != NULL) {
        count++;
        current = current->parent;
    }

    return count;
}


static inline Node* createNode(Node *parent,const System& current) {
	Node *newNode = new Node;
    newNode->parent = parent;
    newNode->system = current;
    newNode->ship = parent->ship;
    newNode->distance = calculate_my_distance(current,parent->system);
    if (current.isNeutron) {
        // cerr << "Found neutron" << endl;
        newNode->ship.isSupercharged = true;
    } else {
        newNode->ship.isSupercharged = false;
    }
    newNode->g = parent->g + 1;
    if (current.isScoopable) {
        // cerr << "Scoopable added to the open list at distance " << current.distanceToDestination << " with current at " << parent->system.distanceToDestination << endl;
        newNode->ship.remainingFuel = newNode->ship.base->tankSize;
    } else {
        double fuelCost = calculateFuelUsed(parent->ship,newNode->distance);
        newNode->ship.remainingFuel -= fuelCost;
        double refuelingTime = calculateRefuelingTime(newNode->ship.base->tankSize,fuelCost);
        // cerr << "Non Scoopable added to the open list with refueling time " << refuelingTime << endl;
        newNode->g += refuelingTime;
    }
    /*
    double cost = 0;
    if (newNode->ship.isSupercharged) {
        cost = newNode->system.distanceToDestination / (averageJumpSize(newNode)*2);
    } else {
        cost = newNode->system.distanceToDestination / averageJumpSize(newNode);
    }
    */
    newNode->h = calculate_cost_to_destination(current,newNode->ship);
    newNode->f = newNode->g + newNode->h;
    return newNode;
}

void debugRoute(Node* node) {
    bool first = true;
    cerr << routeSize(node) << " length route with cost " << node->g << " cost route found. (";
    while (node != NULL) {
        if (!first) {
            cerr << ", ";
        }
        first = false;
        cerr << node->system.id64;
        node = node->parent;
    }
    cerr << ")" << endl;
}

int main() {
    System destination;
    destination.id64 = 20578934;
    destination.x = 25.21875;
    destination.y = -20.90625;
    destination.z = 25899.96875;

    System source;
    source.id64 = 10477373803;
    source.x = 0;
    source.y = 0;
    source.z = 0;
    source.isScoopable = false;
    source.isNeutron = false;
    source.distanceToDestination = calculate_my_distance(source,destination);

    ShipDetails shipDetails;
    shipDetails.fuelMultiplier = 0.012;
    shipDetails.fuelPower = 2.6;
    shipDetails.baseMass = 477.85;
    shipDetails.optimalMass = 2901.6;
    shipDetails.tankSize = 32;
    shipDetails.maxFuelPerJump = 8;
    shipDetails.rangeBoost = 10.5;
    shipDetails.scoopSpeed = 1245.00;
    shipDetails.injectionBoost = true;
    Ship ship;
    ship.base = &shipDetails;
    ship.remainingFuel = 32;
    ship.isSupercharged = false;

    /*
    ship.rangeBoost = 0;
    cerr << "Ship has maximum full range of " << calculateMaxRange(ship) << endl;
    ship.rangeBoost = 10.5;
    cerr << "Ship has maximum full range with guardian booster of " << calculateMaxRange(ship) << endl;
    ship.rangeBoost = 0;
    ship.remainingFuel = 8;
    cerr << "Ship has maximum empty range of " << calculateMaxRange(ship) << endl;
    ship.rangeBoost = 10.5;
    cerr << "Ship has maximum empty range with guardian booster of " << calculateMaxRange(ship) << endl;
    ship.rangeBoost = 0;
    ship.remainingFuel = 4;
    cerr << "Ship has maximum half max range of " << calculateMaxRange(ship) << endl;
    ship.rangeBoost = 10.5;
    cerr << "Ship has maximum half max range with guardian booster of " << calculateMaxRange(ship) << endl;

    ship.remainingFuel = 64;
    double distance = calculateMaxRange(ship);
    cerr << "Ship would use " << calculateFuelUsed(ship,distance) << " tonnes of fuel to travel max range" << endl;
    cerr << "Ship would take " << calculateRefuelingTime(ship,distance) << " seconds to refuel that" << endl;
    distance = 10;
    cerr << "Ship would use " << calculateFuelUsed(ship,distance) << " tonnes of fuel to travel " << distance << "LY" << endl;
    cerr << "Ship would take " << calculateRefuelingTime(ship,distance) << " seconds to refuel that" << endl;

    ship.remainingFuel = 0;
    cerr << "Ship would take " << calculateRefuelingTime(ship) << " seconds to fully refuel from " << ship.remainingFuel << " tonnes left" << endl;
    ship.remainingFuel = 32;
    cerr << "Ship would take " << calculateRefuelingTime(ship) << " seconds to fully refuel from " << ship.remainingFuel << " tonnes left" << endl;
    ship.remainingFuel = 64;
    */

    Node sourceNode;
    sourceNode.parent = NULL;
    sourceNode.system = source;
    sourceNode.ship = ship;
    sourceNode.distance = 0;
    sourceNode.g = 0;
    sourceNode.h = calculate_cost_to_destination(source,ship);
    sourceNode.f = sourceNode.h;

	auto compare = [] (const Node* a, const Node* b) {
        // :cerr << "Comparing " << a->system.id64 << "(" << a->f << ") to " << b->system.id64 << "(" << b->f << ")" << endl;
        return a->f > b->f;
    };

    VisitedNode closedMap;
    closedMap.set_empty_key(0);
    closedMap.set_deleted_key(0xffffffff);

    __gnu_pbds::priority_queue<Node*,decltype(compare),__gnu_pbds::binomial_heap_tag> open(compare);

    open.push(&sourceNode);

    SystemFetcher fetcher;

    unsigned long int maximumSize = 0;

    do {
        Node *current = open.top();
        open.pop();

        VisitedNode::iterator isClosedItr = closedMap.find(current->system.id64);
        if (isClosedItr != closedMap.end()) {
            continue;
        }
        closedMap[current->system.id64] = current;


        string isNeutron = "";
        if (current->system.isNeutron) {
            isNeutron = "neutron ";
        }
        string isScoopable = "";
        if (current->system.isScoopable) {
            isScoopable = "scoopable ";
        }

        unsigned long int currentSize = routeSize(current);

        if (currentSize >= maximumSize) {
            maximumSize = currentSize;
            double distance = calculate_my_distance(current->system,destination);
            cerr << "Route size " << currentSize << " going " << current->distance << "LY to current " << isNeutron << isScoopable << current->system.id64 << "(" << current->system.x << "," << current->system.y << "," << current->system.z << ") at distance " << distance << "LY at " << (source.distanceToDestination-distance)/currentSize << "LY/jump with fuel left " << current->ship.remainingFuel << " tonnes and with cost " << current->g << " and estimated cost " << current->f << " and " << open.size() << " items on the open list and " << closedMap.size() << " items on the closed list" << endl;
            // debugRoute(current);
            //if (open.size() > 0 && open.top()->ship.remainingFuel < open.top()->ship.base->maxFuelPerJump) {
            //    cerr << "Next on list has cost " << open.top()->f << " and distance " <<  calculate_my_distance(open.top()->system,destination) << " and is " << open.top()->system.id64 << endl;
            //}
        }

        double currentRange = calculateMaxRange(current->ship);
        if (current->system.distanceToDestination <= currentRange) {
            cerr << current->system.distanceToDestination << " " << currentRange << endl;
            debugRoute(current);
            break;
        }
        // cerr << current->system.x << "," << current->system.y << "," << current->system.z << endl;
        // cerr << destination.x << "," << destination.y << "," << destination.z << endl;
        vector<System> systems = fetcher.lookup(current->system,destination,currentRange);

        vector<Node*> newOpen;

        VisitedNode toDelete;
        toDelete.set_empty_key(0);
        toDelete.set_deleted_key(0xffffffff);

        for (const System& system : systems) {
            unsigned long int id64 = system.id64;
            VisitedNode::iterator isClosedItr = closedMap.find(id64);
            if (isClosedItr != closedMap.end()) {
                // cerr << "Skipping " << id64 << " as it is already on the closed list" << endl;
                continue;
            }
            // cerr << "Continuing with " << id64 << endl;
            Node* newNode = createNode(current,system);

            newOpen.push_back(newNode);
            // cout << "Pushed " << myRow["id64"].as<unsigned long int>() << " was " << myRow["distance_from_source"] << "LY away from source and " << myRow["distance_to_destination"] << "LY away from destination" << endl;
        }
        /*
        if (toDelete.size()) {
            while (open.size()) {
                Node *current = open.top();
                open.pop();
                VisitedNode::iterator isDeletedItr = toDelete.find(current->system.id64);
                if (isDeletedItr == toDelete.end()) {
                    newOpen.push_back(current);
                } else {
                    delete current;
                }
            }
        }
        */
        for (auto node : newOpen) {
            open.push(node);
        }
    } while (open.size() > 0 && (LIMIT == 0 || open.size() < LIMIT));

    while (open.size() > 0) {
        Node *current = open.top();
        open.pop();
        delete current;
    }
}
