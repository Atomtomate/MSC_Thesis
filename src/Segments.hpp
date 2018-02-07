#ifndef SEGMENTS_HPP_
#define SEGMENTS_HPP_

#include "Config.hpp"

#include <string>
#include <set>
#include <limits>

namespace DMFT
{
    template<unsigned char NFLAVORS>
    class Segments
    {
        public:
            typedef std::pair<const RealT, const RealT> SConfig;
            struct compare {
                inline bool operator()(const SConfig& lhs, const SConfig& rhs) const
                {
                    return (lhs.first) < (rhs.first);
                }
            };

            Segments(const Config &conf): conf(&conf), fullLine({})
            {
                for(int i = 0; i < NFLAVORS; i++) fullLine[i] = false;
            }

            /*! Computed correct time ordered list of segments for periodic bc
             *
             *  @param [in] list list of pairs containing start and end points
             *              of segments
             *  @return list of segments with correct ordering under periodic bc
             */
            std::vector<std::pair<RealT,RealT> > timeOrder(std::initializer_list<std::pair<RealT,RealT> > list) const
            {
                std::vector<std::pair<RealT,RealT> > res;
                res.reserve(list.size());
                bool wrapT = false;
                for(auto el : list)
                {
                    if(el.first > conf->beta) wrapT = true;
                    if(wrapT)
                    {
                        if(std::fmod(el.first,conf->beta) >= std::fmod(el.second,conf->beta))
                            LOG(WARNING) << "Too many segments requested, segments wrap around more than one time!";
                        res.emplace_back(std::fmod(el.first,conf->beta), std::fmod(el.second,conf->beta));
                    }
                    else
                    {
                        if(el.second >= conf->beta) wrapT = true;
                        if(el.first >= el.second)
                        {
                            wrapT = true;
                            if(el.second >= conf->beta) LOG(WARNING) << "Segments wrap around more than one time!";
                            res.emplace_back(el.first, el.second + conf->beta);
                        }
                        else
                            res.emplace_back(el.first, el.second);
                    }
                }
                return res;
            }

            int sign(const unsigned int f) const
            {
                if(fullLine[f] || lines[f].size() == 0) return 1;
                auto it = lines[f].rbegin();
                if(it->second > conf->beta)
                    return -1;
                return 1;
            }

            const SConfig getTimeOrdered(unsigned int i, const unsigned int f) const
            {
                i = i % lines[f].size();
                auto it = lines[f].cbegin();
                std::advance(it, i);
                const SConfig res = *it;
                return res; 
            }

            const int getIndex(const RealT ts, const unsigned int f)
            {
                unsigned int index = -1;
                for(auto el : lines[f])
                {
                    index += 1;
                    if(inSegment(ts, el.first, el.second))
                        return index;
                }
                return -1;
            }

            const SConfig getByT(const RealT ts, const unsigned int f)
            {
                return getTimeOrdered(getIndex(ts, f), f);
            }

            bool inSegment(RealT t, const RealT ts, const RealT tf)
            {
                if(( ts <= t && tf > t ) || ( tf > conf->beta && (tf - conf->beta) > t))
                    return true;
                t = std::fmod(t, conf->beta); 
                if(( ts <= t && tf > t ) || ( tf > conf->beta && (tf - conf->beta) > t))
                    return true;
                return false;
            }

            RealT overlap()
            {
                RealT res = 0.;
                for(int f = 0; f < (int)((NFLAVORS+1)/2); f++)
                {
                    for(auto el : lines[f])
                    {
                        res += overlap(el.first, el.second, f);
                    }
                }
                return res;
            }
            
            RealT overlap(RealT ts, RealT tf, unsigned char f)
            {
                RealT res = 0.;
                if(tf < ts) tf += conf->beta;
                for(int fp = 0; fp < NFLAVORS; fp++)
                {
                    if(f == fp) continue;
                    if(fullLine[fp])
                    {
                        res += tf-ts;
                        continue;
                    }
                    for(auto el : lines[fp])
                    {
                        res += std::max(0., std::min(el.second, tf) - std::max(el.first, ts));
                        if(el.second > conf->beta)                              // case: [ ==  el.f --   el.s ==]
                            res += std::max(0., std::min(el.second - conf->beta, tf) - ts);
                        if(tf > conf->beta)
                            res += std::max(0., std::min(el.second, tf - conf->beta)  - el.first );
                    }
                }
                return res;
            }


            /*! returns the amount of empty line after starting point. negative
             * if within segment
             *
             *  @param [in] ts  point on line
             *  @param [in] f   flavor
             *
             *  @return amount of empty time line behind point
             */
            RealT maxl(const RealT ts, const unsigned int f) const
            {
                // border cases
                if(lines[f].size() == 0) return conf->beta; 
                if(fullLine[f]) return -conf->beta;
                auto it0 = lines[f].rbegin();
                if(ts < (it0->second - conf->beta))                             // ts falls within last segment that wraps around
                    return -(it0->second - conf->beta - ts);
                auto it = lines[f].cbegin();
                if(ts < it->first)
                    return it->first-ts;                                        // ts falls before the first segment
                // otherwise:
                while(it != lines[f].cend())
                {
                    auto tsi = it->first;
                    auto tfi = it->second;
                    if(ts == tsi) return tsi - tfi;
                    if(ts > tsi && ts < tfi) return -(tfi - ts);                // ts falls within segment
                    ++it;
                    if(ts > tfi)                                                // in free space behind segment
                    {   
                        if(it == lines[f].cend())                              // ts falls behind last segment
                        {
                           it = lines[f].cbegin();
                           return (it->first + conf->beta - ts);
                        }
                        else if(ts < it->first)                                // ts falls in free space
                            return (it->first - ts);
                    }
                }
                LOG(WARNING) << "Maximum length could not be computed. ts == ts1 in list. This is unlikely. f = " << f << ", ts = " << ts;
                return 0;
            }

            RealT dist_to_next_end(const RealT ts, const unsigned int f)
            {
                if(lines[f].size() == 0 || fullLine[f]) return conf->beta;
                else if(lines[f].size() == 1)
                {
                    auto res = lines[f].cbegin()->second - ts;
                    res = res < 0 ? res + conf->beta : res;
                    return res;
                }
                auto ml = maxl(ts, f);
                int index = 0;
                if(ml < 0)  // already in segment, get indexof next one
                    index = getIndex(ts, f) + 1;
                else        // between egments, forward to next segment and get index
                    index = getIndex(std::nextafter(ts+ml, 2*conf->beta), f);
                RealT res = getTimeOrdered(index, f).second - ts;
                res = res <= 0 ? res + conf->beta : res;
                return res;
            }


            bool insertFullLine(const unsigned int f)
            {
                if(fullLine[f] || lines[f].size() != 0) return false;
                auto res = lines[f].insert(std::make_pair(0., std::nextafter(conf->beta,0)));
                if(res.second > 0)
                {
                    fullLine[f] = true;
                    return true;
                }
                LOG(WARNING) << "Full line insertion failed.";
                return false;
            }

            inline const bool hasFullLine(const unsigned int f) const
            {
                return fullLine[f];
            }

            int insertSegment(const RealT ts, const RealT tf, const unsigned int f)
            {
                if(maxl(ts, f) < (tf-ts) || tf < ts) return -1;                   // no valid sgment
                if(fullLine[f])
                {
                    LOG(WARNING) << "Trying to insert segment on full line.";
                    return -1;
                }
                else if(std::abs(conf->beta - tf + ts) < 2*std::numeric_limits<RealT>::epsilon())
                {
                    LOG(WARNING) << "Trying to insert >almost< full line. Forcing full line insertion!";
                    insertFullLine(f);
                    return 0;
                }
                else if(tf - ts < 2*std::numeric_limits<RealT>::epsilon())
                {
                    LOG(WARNING) << "Trying to insert >almost< zero segment. Rejecting insertion!";
                    return -1;
                }
                auto res = lines[f].insert(std::make_pair(ts, tf));
                if(res.second > 0)
                    return std::distance(lines[f].begin(), res.first);
                return -1;
            }

            SConfig deleteSegment(const unsigned int i, const unsigned char f)
            {
                if(i >= lines[f].size())
                {
                    LOG(WARNING) << "Trying to delete invalid segment.";
                    return std::make_pair(-1.,-1.);
                }
                if(fullLine[f])
                {
                    if(lines[f].size() != 1)
                    {
                        LOG(ERROR) << "Found full line, but more than one segment!";
                        return std::make_pair(-1.,-1.);
                    }
                    fullLine[f] = false;
                    auto it = lines[f].begin();
                    auto c = *it;
                    lines[f].erase(it);
                    return c;
                }
                else
                {
                    auto it = lines[f].begin();
                    std::advance (it,i);
                    auto c = *it;
                    lines[f].erase(it);
                    return c;
                }
            }

            SConfig deleteSegmentT(const RealT ts, const unsigned char f)
            { 
                auto it = lines[f].begin();
                if(fullLine[f])
                    return deleteSegment(0, f);
                while(it != lines[f].end())
                {
                    if(ts >= it->first && ts <= it->second)
                    {
                        auto c = *it;
                        lines[f].erase(it);
                        return c;
                    }
                    it++;
                }
                --it;
                if(ts <= it->second - conf->beta)
                {
                    auto c = *it;
                    lines[f].erase(it);
                    return c;
                }
                return std::make_pair(-1,-1);
            }

            /*! Insert anti segment.
             *
             *  @param [in] ts  start of anti segment (end of first segment)
             *  @param [in] tf  end of anti segment (start of second segment)
             *
             *  @return success flag
             */
            bool insertAntiSegment(RealT ts, RealT tf, const unsigned char f)
            {
                auto l = -maxl(ts, f);
                if((l < (tf-ts)) || (l == 0) || (ts >= tf)) return false;
                if(fullLine[f])
                {
                    deleteSegment(0,f);
                    SConfig res;
                    if(tf <= conf->beta)
                        res = insertSegment(tf, ts+conf->beta, f);
                    else
                        res = insertSegment(tf-conf->beta, ts, f);
                    if(res.first < 0)
                    {
                        LOG(ERROR) << "failed to insert anti segment into full line";
                        return false;
                    }
                    return true;
                }
                SConfig old = deleteSegmentT(ts, f);
                if(old.first > 0)
                {
                    ts = ts > old.first ? ts : ts + conf->beta;
                    if(!insertSegment(old.first, ts, f))                         // roll back on failure
                    {
                        LOG(WARNING) << "Error during insertion of first segment in anti segment insertion. ts: " << old.first << " - tf: " << ts;
                        insertSegment(old.first, old.second, f);
                        return false;
                    }
                    tf = tf > conf->beta ? tf - conf->beta : tf;
                    RealT tf1 = old.second;
                    if(tf1 - tf > conf->beta)
                        tf1 = old.second - conf->beta;
                    if(!insertSegment(tf, tf1, f))                        // roll back on failure
                    {
                        LOG(WARNING) << "Error during insertion of second segment in anti segment insertion ts: " << tf << " - tf: " << tf1;
                        deleteSegmentT(old.first, f);
                        insertSegment(old.first, old.second, f);
                        return false;
                    }
                    return true;
                }
                return false;
            }

            bool deleteAntiSegment(const unsigned int i, const unsigned char f)
            {
                SConfig old2 = deleteSegment(i, f);
                if(old2.first < 0)
                {
                    LOG(WARNING) << "Could not remove second old segment for anti segment removal.";
                    return false;
                }
                // Full line at as removal for n_s = 1
                if(lines[f].size() == 0)
                    return insertFullLine(f);
                // more than one segment present, continue normally
                SConfig old1 = deleteSegment( (i > 0 ? i-1 : lines[f].size()-1), f) ;
                if(old1.first < 0)
                {
                    insertSegment(old2.first, old2.second, f);
                    LOG(WARNING) << "Could not remove first old segment for anti segment removal. ts: " << old2.first << " i: " << i;
                    return false;
                }
                LOG(INFO) << " :: " << old1.first << " ---- " << old2.second;
                RealT tf = old1.first < old2.second ? old2.second : old2.second + conf->beta;
                if(!insertSegment(old1.first, tf, f))
                {
                    LOG(WARNING) << "Could not insert new segment of anti segment removal.  ts: " << old1.first << " - tf: " << tf;
                    insertSegment(old1.first, old1.second, f);
                    insertSegment(old2.first, old2.second, f);
                    return false;
                }
                return true;
            }

            std::string print_segments() const
            {
                std::string tmp[NFLAVORS] = {};
                std::string res = "";
                for(unsigned char fi = 0; fi < NFLAVORS; fi++)
                {
                    if(fullLine[fi])
                    {
                        res += "F" + std::to_string(fi) + " 0 == "  + std::to_string(conf->beta) + "\n";
                        if(lines[fi].size() != 1) LOG(WARNING) << "Full line, but not exactly one segment found!";
                        continue;
                    }
                    tmp[fi] +=  " " + std::to_string(fi) + ": ";
                    for(auto it = lines[fi].begin(); it != lines[fi].end(); it++)
                    {
                        tmp[fi] += " -- " + std::to_string(it->first) + "=="  + std::to_string(it->second);
                    }
                    res += tmp[fi] + "\n";
                }
                return res;
            }
            
        private:
            typedef std::set<SConfig, compare> SConfigL;                                  // Time orderer list (set)
            Config const * conf;
            std::array<bool, NFLAVORS> fullLine;
            SConfigL lines[NFLAVORS];
    }; 
}
#endif
