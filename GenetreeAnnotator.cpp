#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <unordered_set>
#include <unordered_map>

using namespace std;

void display(vector<pair<string, int> > tr){
	for (auto e: tr) cerr << e.second << "\t" << e.first << endl;
	cerr << endl;
}

vector<string> sample(vector<pair<string, int> > tr, int level){
	vector<string> result;
	vector<pair<string, int> > stack;
	for (auto &e: tr){
		if (e.second <= level){
			for (int i = 0, n = 1 << (level - e.second); i < n; i++) result.push_back(e.first);
			continue;
		}
		stack.push_back(e);
		while (stack.size() >= 2 && stack[stack.size() - 2].second == stack.back().second){
			if (rand() & 1) stack[stack.size() - 2].first = stack.back().first;
			stack[stack.size() - 2].second--;
			stack.pop_back();
		}
		if (stack.back().second == level){
			result.push_back(stack.back().first);
			stack.pop_back();
		}
	}
	return result;
}

class DynamicBitset{
	int size = 0;
	vector<uint64_t> vec;
	
public:
	DynamicBitset(){}
	DynamicBitset(int sz): size(sz), vec((sz + 63) / 64){}
	
	void set(int i){
		if (i >= size){
			size = i + 1;
			if ((size + 63) / 64 > vec.size()){
				vec.resize((size + 63) / 64);
			}
		}
		vec[i / 64] |= (1LL << (i % 64));
	}
	
	DynamicBitset operator|(const DynamicBitset &b) const{
		if (size < b.size) return b | *this;
		DynamicBitset res(size);
		for (int i = 0; i < b.vec.size(); i++){
			res.vec[i] = vec[i] | b.vec[i];
		}
		for (int i = b.vec.size(); i < vec.size(); i++){
			res.vec[i] = vec[i];
		}
		return res;
	}
	
	DynamicBitset operator&(const DynamicBitset &b) const{
		if (size < b.size) return b & *this;
		DynamicBitset res(b.size);
		for (int i = 0; i < b.vec.size(); i++){
			res.vec[i] = vec[i] & b.vec[i];
		}
		return res;
	}
	
	DynamicBitset operator^(const DynamicBitset &b) const{
		if (size < b.size) return b ^ *this;
		DynamicBitset res(size);
		for (int i = 0; i < b.vec.size(); i++){
			res.vec[i] = vec[i] ^ b.vec[i];
		}
		for (int i = b.vec.size(); i < vec.size(); i++){
			res.vec[i] = vec[i];
		}
		return res;
	}
	
	DynamicBitset operator-(const DynamicBitset &b) const{
		DynamicBitset res(size);
		for (int i = 0; i < vec.size(); i++){
			if (i < b.vec.size()) res.vec[i] = vec[i] & ~b.vec[i];
			else res.vec[i] = vec[i];
		}
		return res;
	}
	
	bool operator==(const DynamicBitset &b) const{
		if (size < b.size) return b == *this;
		for (int i = 0; i < b.vec.size(); i++){
			if (vec[i] != b.vec[i]) return false;
		}
		for (int i = b.vec.size(); i < vec.size(); i++){
			if (vec[i] != 0) return false;
		}
		return true;
	}
	
	bool operator!=(const DynamicBitset &b) const{
		return !(*this == b);
	}
	
	bool isDisjointTo(const DynamicBitset &b) const{
		if (size < b.size) return b.isDisjointTo(*this);
		for (int i = 0; i < b.vec.size(); i++){
			if ((vec[i] & b.vec[i]) != 0) return false;
		}
		return true;
	}
	
	vector<int> setBits() const{
		vector<int> res;
		for (int i = 0; i < vec.size(); i++){
			for (int j = 0; j < 64; j++){
				if (vec[i] & (1LL << j)) res.push_back(64 * i + j);
			}
		}
		return res;
	}
};

struct ClusterHash{
	uint64_t h[2] = {};
	
	ClusterHash(){}
	
	ClusterHash(uint64_t a, uint64_t b){
		h[0] = a;
		h[1] = b;
	}
	
	ClusterHash operator + (const ClusterHash &c) const{
		return ClusterHash(h[0] + c.h[0], h[1] + c.h[1]);
	}
	
	ClusterHash operator - (const ClusterHash &c) const{
		return ClusterHash(h[0] - c.h[0], h[1] - c.h[1]);
	}
	
	bool operator == (const ClusterHash &c) const{
		return h[0] == c.h[0] && h[1] == c.h[1];
	}
	
	bool operator < (const ClusterHash &c) const{
		return h[0] < c.h[0] || (h[0] == c.h[0] && h[1] < c.h[1]);
	}
};

struct ClusterHasher{
	size_t operator()(const ClusterHash &c) const{
		return c.h[0];
	}
};

struct ClusterHashLinkedListNode{
	ClusterHash c;
	ClusterHashLinkedListNode *next = nullptr;
	~ClusterHashLinkedListNode(){
		if (next != nullptr) delete next;
	}
	
	bool operator == (const ClusterHashLinkedListNode &n) const{
		if (this == nullptr && &n == nullptr) return true;
		if (this == nullptr || &n == nullptr) return false;
		return c == n.c && *next == *n.next;
	}
};

struct PartitionHash{
	ClusterHash y, x0, x1;
	shared_ptr<ClusterHashLinkedListNode> xs;
	
	PartitionHash(){}
	
	PartitionHash(const ClusterHash &y, const ClusterHash &x0, const ClusterHash &x1): y(y), x0(min(x0, x1)), x1(max(x0, x1)){}
	
	PartitionHash(const ClusterHash &y, vector<ClusterHash> &xlist): y(y){
		sort(xlist.begin(), xlist.end());
		x0 = xlist[0];
		x1 = xlist[1];
		if (xlist.size() > 2){
			ClusterHashLinkedListNode *p = new ClusterHashLinkedListNode;
			p->c = xlist[2];
			xs.reset(p);
			for (int i = 3; i < xlist.size(); i++){
				p->next = new ClusterHashLinkedListNode;
				p = p->next;
				p->c = xlist[i];
			}
		}
	}
	
	bool operator == (const PartitionHash &p) const{
		return y == p.y && x0 == p.x0 && x1 == p.x1 && *xs == *(p.xs);
	}
	
	size_t hash() const{
		uint64_t res = y.h[0] + x0.h[0] * x0.h[1] + x1.h[0] * x1.h[1];
		for (ClusterHashLinkedListNode* p = xs.get(); p != nullptr; p = p->next){
			res += p->c.h[0] * p->c.h[1]; 
		}
		return res;
	}
};

struct PartitionHasher{
	size_t operator()(const PartitionHash &p) const{
		uint64_t res = p.y.h[0] + p.x0.h[0] * p.x0.h[1] + p.x1.h[0] * p.x1.h[1];
		for (ClusterHashLinkedListNode* q = p.xs.get(); q != nullptr; q = q->next){
			res += q->c.h[0] * q->c.h[1]; 
		}
		return res;
	}
};

class GenetreeAnnotator{
public:
	struct Polytree{
		vector<ClusterHash> singletons;
		unordered_set<ClusterHash, ClusterHasher> clusters;
		unordered_set<PartitionHash, PartitionHasher> additions;
		unordered_map<PartitionHash, int, PartitionHasher> partitions;
		
		void write(string file, const vector<string> &id2name) const{
			ofstream fout(file);
			fout << singletons.size() << endl;
			for (const string &s: id2name) fout << s << endl;
			unordered_map<ClusterHash, int, ClusterHasher> clusterId;
			int idCounter = 0;
			for (auto &c: singletons) clusterId[c] = idCounter++;
			for (auto &c: clusters){
				if (clusterId.count(c) == 0) clusterId[c] = idCounter++;
			}
			fout << clusterId.size() << "\t" << additions.size() << endl;
			for (auto &p: additions) fout << clusterId[p.y] << "\t" << clusterId[p.x0] << "\t" << clusterId[p.x1] << endl;
			fout << partitions.size() << endl;
			for (auto &p: partitions) fout << clusterId[p.first.y] << "\t" << clusterId[p.first.x0] << "\t" << clusterId[p.first.x1] << "\t" << p.second << endl;
		}
	};
private:
	struct Node{
		int leftChildId = -1, rightChildId = -1;
		DynamicBitset label;
		bool isDuplication = false;
		bool isLeaf = false;
		int score = -1;
		int leafId = -1;
	};
	
	vector<Node> node;
	vector<string> id2name;
	unordered_map<string, int> name2id;
	
	tuple<int, int, int> createSubtree(const unordered_map<long long, string> &leafname, const unordered_map<long long, pair<long long, long long> > &children, const long long cur, vector<int> &rootId){
		if (children.count(cur) == 0){
			int curId = node.size();
			node.emplace_back();
			node[curId].isLeaf = true;
			node[curId].score = 0;
			if (name2id.count(leafname.at(cur)) == 0){
				name2id[leafname.at(cur)] = id2name.size();
				id2name.push_back(leafname.at(cur));
			}
			node[curId].leafId = name2id[leafname.at(cur)];
			node[curId].label.set(name2id[leafname.at(cur)]);
			return make_tuple(curId, -1, -1);
		}
		
		tuple<int, int, int> left = createSubtree(leafname, children, children.at(cur).first, rootId), right = createSubtree(leafname, children, children.at(cur).second, rootId);
		int cur0 = node.size();
		node.emplace_back();
		int cur1 = node.size();
		node.emplace_back();
		int cur2 = node.size();
		node.emplace_back();
		node[cur0].leftChildId = get<0>(left);
		node[cur0].rightChildId = get<0>(right);
		
		node[cur1].leftChildId = get<0>(left);
		node[cur2].leftChildId = get<0>(right);
		if (get<1>(left) != -1) node[get<1>(left)].rightChildId = cur2;
		if (get<2>(left) != -1) node[get<2>(left)].rightChildId = cur2;
		if (get<1>(right) != -1) node[get<1>(right)].rightChildId = cur1;
		if (get<2>(right) != -1) node[get<2>(right)].rightChildId = cur1;
		
		int root1 = node.size();
		node.emplace_back();
		node[root1].leftChildId = get<0>(left);
		node[root1].rightChildId = cur2;
		rootId.push_back(root1);
		int root2 = node.size();
		node.emplace_back();
		node[root2].leftChildId = get<0>(right);
		node[root2].rightChildId = cur1;
		rootId.push_back(root2);
		
		return make_tuple(cur0, cur1, cur2);
	}
	
	int scoreSubtree(int cur){
		if (node[cur].score != -1) return node[cur].score;
		node[cur].score = scoreSubtree(node[cur].leftChildId) + scoreSubtree(node[cur].rightChildId);
		node[cur].label = node[node[cur].leftChildId].label | node[node[cur].rightChildId].label;
		if (!node[node[cur].leftChildId].label.isDisjointTo(node[node[cur].rightChildId].label)){
			node[cur].score++;
			node[cur].isDuplication = true;
			if (node[cur].label != node[node[cur].leftChildId].label) node[cur].score++;
			if (node[cur].label != node[node[cur].rightChildId].label) node[cur].score++;
		}
		return node[cur].score;
	}
	
	static vector<pair<string, int> > merge(const vector<pair<string, int> > &a, const vector<pair<string, int> > &b){
		vector<pair<string, int> > c;
		int p = 0, q = 0, i = -1, j = -1;
		while (i != a.size() - 1 || j != b.size() - 1){
			if (p == 0){
				i++;
				j++;
				q = max(a[i].second, b[j].second);
				p = (1 << (q - a[i].second)) - (1 << (q - b[j].second));
				c.push_back({"(" + a[i].first + "," + b[j].first + ")", q});
			}
			else if (p < 0){
				i++;
				int nq = max(a[i].second, q);
				int np = (p << (nq - q)) + (1 << (nq - a[i].second));
				p = np;
				q = nq;
				c.push_back({"(" + a[i].first + "," + b[j].first + ")", a[i].second});
			}
			else {
				j++;
				int nq = max(b[j].second, q);
				int np = (p << (nq - q)) - (1 << (nq - b[j].second));
				p = np;
				q = nq;
				c.push_back({"(" + a[i].first + "," + b[j].first + ")", b[j].second});
			}
		}
		return c;
	}
	
	static bool eq(const vector<pair<string, int> > &a, const vector<pair<string, int> > &b){
		if (a.size() != b.size()) return false;
		for (int i = 0; i < a.size(); i++){
			if (a[i].first != b[i].first || a[i].second != b[i].second) return false;
		}
		return true;
	}
	
	pair<vector<pair<string, int> >, vector<pair<string, int> > > breakSubtree(int cur) const{
		if (node[cur].isLeaf){
			string name = id2name[node[cur].leafId];
			return {vector<pair<string, int> >(1, {name, 0}), vector<pair<string, int> >(1, {name, 0})};
		}
		auto left = breakSubtree(node[cur].leftChildId), right = breakSubtree(node[cur].rightChildId);
		if (node[cur].isDuplication){
			vector<pair<string, int> > first, second;
			if (eq(left.first, left.second)){
				first = left.first;
			}
			else {
				for (auto &e: left.first) first.push_back({e.first, e.second + 1});
				for (auto &e: left.second) first.push_back({e.first, e.second + 1});
			}
			if (eq(right.first, right.second)){
				second = right.first;
			}
			else {
				for (auto &e: right.first) second.push_back({e.first, e.second + 1});
				for (auto &e: right.second) second.push_back({e.first, e.second + 1});
			}
			return {first, second};
		}
		else {
			return {merge(left.first, right.first), merge(left.second, right.second)};
		}
	}

	ClusterHash buildPolytreePre(int cur, unordered_map<int, ClusterHash> &nodeCluster, Polytree& pt, bool isRoot = false) const{
		if (node[cur].isLeaf) return nodeCluster[cur] = pt.singletons[node[cur].leafId];
		
		int left = node[cur].leftChildId, right = node[cur].rightChildId;
		if (node[cur].isDuplication){
			ClusterHash lc = buildPolytreePre(left, nodeCluster, pt, isRoot), rc = buildPolytreePre(right, nodeCluster, pt, isRoot);
			cerr << "Duplication Called\n";
			if (node[left].label == node[cur].label) return nodeCluster[cur] = lc;
			if (node[right].label == node[cur].label) return nodeCluster[cur] = rc;
			ClusterHash c = lc;
			for (int i: (node[cur].label - node[left].label).setBits()){
				ClusterHash nc = c + pt.singletons[i];
				pt.clusters.insert(nc);
				pt.additions.insert(PartitionHash(nc, c, pt.singletons[i]));
				c = nc;
			}
			return nodeCluster[cur] = c;
		}
		else {
			ClusterHash lc = buildPolytreePre(left, nodeCluster, pt), rc = buildPolytreePre(right, nodeCluster, pt);
			ClusterHash c = lc + rc;
			if (!isRoot) pt.clusters.insert(c);
			if (!isRoot) pt.additions.insert(PartitionHash(c, lc, rc));
			return nodeCluster[cur] = c;
		}
	}
	
	void buildPolytreePost(int cur, unordered_map<int, ClusterHash> &nodeCluster, Polytree& pt, const ClusterHash &s, bool isRoot = false) const{
		if (node[cur].isLeaf) return;
		
		int left = node[cur].leftChildId, right = node[cur].rightChildId;
		if (node[cur].isDuplication){
			buildPolytreePost(left, nodeCluster, pt, s, isRoot);
			buildPolytreePost(right, nodeCluster, pt, s, isRoot);
		}
		else {
			ClusterHash lc = nodeCluster[left], rc = nodeCluster[right];
			if (!isRoot) pt.partitions[PartitionHash(s, lc, rc)]++;
			ClusterHash ls = s + rc;
			pt.clusters.insert(ls);
			if (!isRoot) pt.additions.insert(PartitionHash(ls, s, rc));
			buildPolytreePost(left, nodeCluster, pt, ls);
			ClusterHash rs = s + lc;
			pt.clusters.insert(rs);
			if (!isRoot) pt.additions.insert(PartitionHash(rs, s, lc));
			buildPolytreePost(right, nodeCluster, pt, rs);
		}
	}
		
public:
	const vector<string> &leafnames() const{
		return id2name;
	}
	
	int annotateTree(const unordered_map<long long, string> &leafname, const unordered_map<long long, pair<long long, long long> > &children, const long long root){
		vector<int> rootId;
		tuple<int, int, int> left = createSubtree(leafname, children, children.at(root).first, rootId), right = createSubtree(leafname, children, children.at(root).second, rootId);
		if (get<1>(left) != -1) node[get<1>(left)].rightChildId = get<0>(right);
		if (get<2>(left) != -1) node[get<2>(left)].rightChildId = get<0>(right);
		if (get<1>(right) != -1) node[get<1>(right)].rightChildId = get<0>(left);
		if (get<2>(right) != -1) node[get<2>(right)].rightChildId = get<0>(left);
		
		int curId = node.size();
		node.emplace_back();
		node[curId].leftChildId = get<0>(left);
		node[curId].rightChildId = get<0>(right);
		rootId.push_back(curId);
		
		int bestscore = 999999, bestroot = -1;
		for (int root: rootId){
			int score = scoreSubtree(root);
			if (score < bestscore){
				bestscore = score;
				bestroot = root;
			}
		}
		return bestroot;
	}
	
	vector<pair<string, int> > breakGenetree(int root) const{
		pair<vector<pair<string, int> >, vector<pair<string, int> > > result = breakSubtree(root);
		vector<pair<string, int> > ret;
		if (eq(result.first, result.second)){
			for (auto &e: result.first) ret.push_back({e.first + ";", e.second});
		}
		else {
			for (auto &e: result.first) ret.push_back({e.first + ";", e.second + 1});
			for (auto &e: result.second) ret.push_back({e.first + ";", e.second + 1});
		}
		return ret;
	}
	
	void buildPolytree(int root, Polytree& pt) const{
		while (pt.singletons.size() < id2name.size()){
			ClusterHash c((((uint64_t) rand()) << 60) ^ (((uint64_t) rand()) << 45) ^ (((uint64_t) rand()) << 30) ^ (((uint64_t) rand()) << 15) ^ ((uint64_t) rand()),
				(((uint64_t) rand()) << 60) ^ (((uint64_t) rand()) << 45) ^ (((uint64_t) rand()) << 30) ^ (((uint64_t) rand()) << 15) ^ ((uint64_t) rand()));
			pt.clusters.insert(c);
			pt.singletons.push_back(c);
		}
		unordered_map<int, ClusterHash> nodeCluster;
		buildPolytreePre(root, nodeCluster, pt, true);
		ClusterHash s;	
		buildPolytreePost(root, nodeCluster, pt, s, true);
	}
};

int main(int argc, char** argv) {
	ifstream fin("genetrees.txt");
	ofstream fout("breaked_trees.tre");
	ofstream fx("cluster_trees.tre");
	int k;
	fin >> k;
	GenetreeAnnotator ga;
	GenetreeAnnotator::Polytree pt;
	for (int i = 0; i < k; i++){
		int n, nc;
		long long root;
		fin >> n >> nc >> root;
		unordered_map<long long, string> leafname;
		unordered_map<long long, pair<long long, long long> > children;
		for (int j = 0; j < n; j++){
			long long node;
			string name;
			fin >> name >> node;
			leafname[node] = name;
		}
		for (int j = 0; j < nc; j++){
			long long p, c1, c2;
			fin >> p >> c1 >> c2;
			children[p] = {c1, c2};
		}
		int iroot = ga.annotateTree(leafname, children, root);
		ga.buildPolytree(iroot, pt);
		
		auto trees = ga.breakGenetree(iroot);
		for (auto e: trees) if (count(e.first.begin(), e.first.end(), ',') >= 3) fx << e.first << endl;
		for (string s: sample(trees, 4)) if (count(s.begin(), s.end(), ',') >= 3) fout << s << endl;
	}
	pt.write("polytree.txt", ga.leafnames());
	return 0;
}
