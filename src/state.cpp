#include "state.hpp"
#include <fstream>
#include <type_traits>

static_assert(std::is_trivially_copyable<fac_entry>::value, "fac_entry is written as raw bytes");

static const uint64_t MAGIC_STATE = 0x31455441545346ULL; // "FSTATE1"

template <class T> static void put(std::ostream& o, const T& x) { o.write(reinterpret_cast<const char*>(&x), sizeof x); }
template <class T> static bool get(std::istream& i, T& x) { i.read(reinterpret_cast<char*>(&x), sizeof x); return (bool)i; }

static void put_str(std::ostream& o, const std::string& s) { put<uint64_t>(o, s.size()); o.write(s.data(), (std::streamsize)s.size()); }
static bool get_str(std::istream& i, std::string& s)
{
	uint64_t len;
	if (!get(i, len) || len > (uint64_t(1) << 32)) return false;
	s.resize(len);
	if (len) i.read(&s[0], (std::streamsize)len);
	return (bool)i;
}

template <class T> static void put_vec(std::ostream& o, const vector<T>& v)
{
	put<uint64_t>(o, v.size());
	if (!v.empty()) o.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(T)));
}
template <class T> static bool get_vec(std::istream& i, vector<T>& v)
{
	uint64_t len;
	if (!get(i, len) || len > (uint64_t(1) << 32)) return false;
	v.resize(len);
	if (len) i.read(reinterpret_cast<char*>(v.data()), (std::streamsize)(len * sizeof(T)));
	return (bool)i;
}

static void put_rel(std::ostream& o, const nfs_relation& r)
{
	put(o, r.a); put(o, r.b);
	put_vec(o, r.cols); put_vec(o, r.rfac); put_vec(o, r.afac);
	put(o, r.lp_r); put(o, r.lp_a); put(o, r.lp_a_root);
	put_vec(o, r.members);
}

static bool get_rel(std::istream& i, nfs_relation& r)
{
	return get(i, r.a) && get(i, r.b) && get_vec(i, r.cols) && get_vec(i, r.rfac) && get_vec(i, r.afac)
	    && get(i, r.lp_r) && get(i, r.lp_a) && get(i, r.lp_a_root) && get_vec(i, r.members);
}

void save_state(const std::string& path, const std::string& n, const std::string& m, const std::string& f_text, uint64_t bound,
                const vector<nfs_relation>& rels, const vector<nfs_relation>& partials, uint64_t mat_rows)
{
	std::ofstream o(path, std::ios::binary);
	if (!o) { cerr << "cannot write " << path << '\n'; return; }
	put<uint64_t>(o, MAGIC_STATE);
	put_str(o, n); put_str(o, m); put_str(o, f_text);
	put<uint64_t>(o, bound);
	put<uint64_t>(o, mat_rows);
	put<uint64_t>(o, rels.size());
	for (const auto& r : rels) put_rel(o, r);
	put<uint64_t>(o, partials.size());
	for (const auto& r : partials) put_rel(o, r);
	o.flush();
	if (!o) cerr << "error while writing " << path << '\n';
	else cerr << "saved " << rels.size() << " full + " << partials.size() << " partial relations to " << path << '\n';
}

bool load_state(const std::string& path, gnfs_state& S)
{
	std::ifstream in(path, std::ios::binary);
	uint64_t magic = 0, nr = 0, np = 0;
	if (!in || !get(in, magic) || magic != MAGIC_STATE) return false;
	S = gnfs_state();
	if (!get_str(in, S.n) || !get_str(in, S.m) || !get_str(in, S.f_text) || !get(in, S.bound) || !get(in, S.mat_rows)) return false;
	if (!get(in, nr) || nr > (uint64_t(1) << 32)) return false;
	S.rels.resize(nr);
	for (auto& r : S.rels) if (!get_rel(in, r)) return false;
	if (!get(in, np) || np > (uint64_t(1) << 32)) return false;
	S.partials.resize(np);
	for (auto& r : S.partials) if (!get_rel(in, r)) return false;
	return true;
}
