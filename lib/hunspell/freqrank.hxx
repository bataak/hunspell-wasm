/* Frequency-based reranking of Hunspell suggestions.
 *
 * Hunspell orders suggestions by the stage that produced them (REP, MAP,
 * swaps, ...) and, inside a stage, by the order of the rules in the .aff
 * file.  It has no notion of how likely a candidate is.  FreqRank adds an
 * optional final stage that reorders the suggestion list by a noisy-channel
 * style score
 *
 *     score(c | x) = d(c, x) - lambda * log P(c)
 *
 * where d is an optimal-string-alignment distance in which a vowel/vowel
 * substitution costs 0.5, a consonant/consonant substitution 1.5, any other
 * substitution, insertion, deletion or adjacent transposition 1, and
 * P(c) = (f(c) + 1) / (N + V) is the add-one estimate from a word-form
 * frequency list.  Lower scores come first; ties keep Hunspell's order.
 * If none of the candidates occurs in the frequency list the list is left
 * as it is: without corpus evidence the distance alone ranks worse than
 * Hunspell's stage order.
 *
 * The frequency file is plain UTF-8 text:
 *     # comment
 *     VOWELS аэоөуүиыеёюя      (characters treated as vowels)
 *     LAMBDA 0.5               (optional, default 0.5)
 *     word<TAB>count           (one surface form per line, lowercase)
 *
 * Nothing changes unless a frequency file is loaded.
 */
#ifndef FREQRANK_HXX_
#define FREQRANK_HXX_

#include <cstdint>
#include <string>
#include <vector>

#include "w_char.hxx"

class FreqRank {
 public:
  // Load a frequency file; returns the number of entries (0 on failure).
  size_t load(const char* path);
  bool loaded() const { return !offsets.empty(); }
  // Reorder slst for the misspelled word (UTF-8 dictionaries only).
  void rank(const std::string& word, std::vector<std::string>& slst, int langnum) const;

 private:
  std::string blob;                  // sorted words, '\0'-terminated
  std::vector<uint32_t> offsets;     // start of each word in blob
  std::vector<uint32_t> counts;      // frequency of each word
  std::vector<unsigned short> vowels;  // sorted UTF-16 code units
  double lambda = 0.5;
  double logden = 0.0;               // log(N + V)

  uint32_t lookup(const std::string& w) const;
  bool is_vowel(unsigned short c) const;
  double distance(const std::vector<unsigned short>& a,
                  const std::vector<unsigned short>& b,
                  const std::vector<unsigned short>& al,
                  const std::vector<unsigned short>& bl) const;
};

#endif
