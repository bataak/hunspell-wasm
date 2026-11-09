/* Frequency-based reranking of Hunspell suggestions; see freqrank.hxx. */
#include "freqrank.hxx"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <utility>

#include "csutil.hxx"

namespace {

std::vector<unsigned short> to_u16(const std::string& s) {
  std::vector<w_char> w;
  u8_u16(w, s);
  std::vector<unsigned short> r;
  r.reserve(w.size());
  for (const auto& c : w)
    r.push_back((unsigned short)c);
  return r;
}

std::string trim(const std::string& s) {
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos)
    return std::string();
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

}  // namespace

size_t FreqRank::load(const char* path) {
  FILE* f = fopen(path, "rb");
  if (!f)
    return 0;
  std::vector<std::pair<std::string, uint32_t> > entries;
  std::vector<unsigned short> vw;
  double lam = 0.5;
  std::string line;
  char buf[4096];
  bool ok = true;
  while (fgets(buf, sizeof(buf), f)) {
    line.assign(buf);
    // join lines longer than the buffer
    while (!line.empty() && line.back() != '\n' && fgets(buf, sizeof(buf), f))
      line.append(buf);
    line = trim(line);
    if (line.empty() || line[0] == '#')
      continue;
    size_t tab = line.find('\t');
    if (tab == std::string::npos) {
      if (line.compare(0, 7, "VOWELS ") == 0) {
        vw = to_u16(trim(line.substr(7)));
      } else if (line.compare(0, 7, "LAMBDA ") == 0) {
        lam = atof(line.c_str() + 7);
      }
      continue;
    }
    std::string w = line.substr(0, tab);
    long n = atol(line.c_str() + tab + 1);
    if (w.empty() || n < 0) {
      ok = false;
      break;
    }
    entries.emplace_back(w, (uint32_t)n);
  }
  fclose(f);
  if (!ok || entries.empty())
    return 0;

  std::sort(entries.begin(), entries.end());
  // merge duplicates
  size_t k = 0;
  for (size_t i = 0; i < entries.size(); ++i) {
    if (k && entries[k - 1].first == entries[i].first)
      entries[k - 1].second += entries[i].second;
    else
      entries[k++] = entries[i];
  }
  entries.resize(k);

  blob.clear();
  offsets.clear();
  counts.clear();
  offsets.reserve(k);
  counts.reserve(k);
  double total = 0.0;
  for (const auto& e : entries) {
    offsets.push_back((uint32_t)blob.size());
    blob.append(e.first);
    blob.push_back('\0');
    counts.push_back(e.second);
    total += e.second;
  }
  std::sort(vw.begin(), vw.end());
  vowels.swap(vw);
  lambda = lam;
  logden = std::log(total + (double)k);
  return k;
}

uint32_t FreqRank::lookup(const std::string& w) const {
  size_t lo = 0, hi = offsets.size();
  while (lo < hi) {
    size_t mid = (lo + hi) / 2;
    int c = strcmp(blob.c_str() + offsets[mid], w.c_str());
    if (c == 0)
      return counts[mid];
    if (c < 0)
      lo = mid + 1;
    else
      hi = mid;
  }
  return 0;
}

bool FreqRank::is_vowel(unsigned short c) const {
  return std::binary_search(vowels.begin(), vowels.end(), c);
}

// Weighted optimal string alignment distance (see freqrank.hxx).
// Vowel classes are tested on lowercased characters; equality is exact.
double FreqRank::distance(const std::vector<unsigned short>& a,
                          const std::vector<unsigned short>& b,
                          const std::vector<unsigned short>& al,
                          const std::vector<unsigned short>& bl) const {
  const size_t n = a.size(), m = b.size();
  std::vector<double> T((n + 1) * (m + 1));
  auto at = [m, &T](size_t i, size_t j) -> double& { return T[i * (m + 1) + j]; };
  for (size_t i = 0; i <= n; ++i)
    at(i, 0) = (double)i;
  for (size_t j = 0; j <= m; ++j)
    at(0, j) = (double)j;
  for (size_t i = 1; i <= n; ++i) {
    for (size_t j = 1; j <= m; ++j) {
      double sub = 0.0;
      if (a[i - 1] != b[j - 1]) {
        bool va = is_vowel(al[i - 1]), vb = is_vowel(bl[j - 1]);
        sub = (va && vb) ? 0.5 : ((!va && !vb) ? 1.5 : 1.0);
      }
      double v = std::min(std::min(at(i - 1, j) + 1.0, at(i, j - 1) + 1.0),
                          at(i - 1, j - 1) + sub);
      if (i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1])
        v = std::min(v, at(i - 2, j - 2) + 1.0);
      at(i, j) = v;
    }
  }
  return at(n, m);
}

namespace {

std::string lowercase(const std::string& s, int langnum) {
  std::vector<w_char> w;
  u8_u16(w, s);
  for (auto& ch : w) {
    unsigned short u = unicodetolower((unsigned short)ch, langnum);
    ch.h = (unsigned char)(u >> 8);
    ch.l = (unsigned char)(u & 0xff);
  }
  std::string r;
  u16_u8(r, w);
  return r;
}

}  // namespace

void FreqRank::rank(const std::string& word, std::vector<std::string>& slst,
                    int langnum, const std::vector<replentry>* reptable) const {
  if (!loaded() || slst.size() < 2)
    return;
  // explicit whole-word REP corrections first, each group ranked on its own
  std::vector<std::string> pinned, rest;
  if (reptable) {
    std::string lw = lowercase(word, langnum);
    std::vector<std::string> outs;
    for (const auto& r : *reptable) {
      if (r.pattern != lw)
        continue;
      for (const auto& o : r.outstrings)
        if (!o.empty())
          outs.push_back(trim(o));
    }
    if (!outs.empty()) {
      for (const auto& c : slst) {
        if (std::find(outs.begin(), outs.end(), lowercase(c, langnum)) != outs.end())
          pinned.push_back(c);
        else
          rest.push_back(c);
      }
    }
  }
  if (pinned.empty()) {
    sort_by_score(word, slst, langnum);
    return;
  }
  sort_by_score(word, pinned, langnum);
  sort_by_score(word, rest, langnum);
  pinned.insert(pinned.end(), rest.begin(), rest.end());
  slst.swap(pinned);
}

void FreqRank::sort_by_score(const std::string& word,
                             std::vector<std::string>& slst,
                             int langnum) const {
  if (slst.size() < 2)
    return;
  // word-form counts first: with no corpus evidence for any candidate,
  // distance alone ranks worse than Hunspell's own order, so keep it
  std::vector<uint32_t> cnt(slst.size());
  bool evidence = false;
  for (size_t i = 0; i < slst.size(); ++i) {
    std::vector<w_char> lw;
    u8_u16(lw, slst[i]);
    for (auto& ch : lw) {
      unsigned short u = unicodetolower((unsigned short)ch, langnum);
      ch.h = (unsigned char)(u >> 8);
      ch.l = (unsigned char)(u & 0xff);
    }
    std::string lower;
    u16_u8(lower, lw);
    cnt[i] = lookup(lower);
    if (cnt[i])
      evidence = true;
  }
  if (!evidence)
    return;

  std::vector<unsigned short> x = to_u16(word);
  std::vector<std::pair<double, size_t> > score;
  score.reserve(slst.size());
  for (size_t i = 0; i < slst.size(); ++i) {
    std::vector<unsigned short> c = to_u16(slst[i]);
    // vowel test uses lowercase; keep exact characters for equality
    std::vector<unsigned short> cl(c), xl(x);
    for (auto& ch : cl) ch = unicodetolower(ch, langnum);
    for (auto& ch : xl) ch = unicodetolower(ch, langnum);
    double d = distance(c, x, cl, xl);
    double lp = std::log((double)cnt[i] + 1.0) - logden;
    score.emplace_back(d - lambda * lp, i);
  }
  std::stable_sort(score.begin(), score.end(),
                   [](const std::pair<double, size_t>& p,
                      const std::pair<double, size_t>& q) { return p.first < q.first; });
  std::vector<std::string> out;
  out.reserve(slst.size());
  for (const auto& s : score)
    out.push_back(slst[s.second]);
  slst.swap(out);
}
