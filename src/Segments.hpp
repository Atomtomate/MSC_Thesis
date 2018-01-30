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

            Segments(const Config &conf): conf(conf) {}


            const SConfig getTimeOrdered(unsigned int i, const unsigned char f) const
            {
                i = i % lines[f].size();
                auto it = lines[f].cbegin();
                std::advance(it, i);
                const SConfig res = *it;
                return res; 
            }

            bool inSegment(const RealT t, const RealT ts, const RealT tf)
            {
                if(( ts < t && tf > t ) || \
                    ( tf > conf.beta &&  tf - conf.beta > t))
                    return true;
                else
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
                if(tf < ts) tf += conf.beta;
                for(int fp = 0; fp < NFLAVORS; fp++)
                {
                    if(f == fp) continue;
                    for(auto el : lines[fp])
                    {
                        if(el.second > conf.beta && tf < conf.beta){
                            res += std::max(0., std::min(el.second, tf + conf.beta) - std::max(el.first, ts + conf.beta));
                        }else if(tf > conf.beta && el.second < conf.beta){
                            res += std::max(0., std::min(el.second + conf.beta, tf)  - std::max(el.first + conf.beta, ts));
                        }else{
                            res += std::max(0., std::min(el.second, tf) - std::max(el.first, ts));
                        }
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
            RealT maxl(const RealT ts, const unsigned char f) const
            {
                // border cases
                if(lines[f].size() == 0) return conf.beta; 
                auto it0 = lines[f].rbegin();
                if(ts < (it0->second - conf.beta))                               // ts falls within last segment that wraps around
                    return -(it0->second - conf.beta - ts);
                auto it = lines[f].cbegin();
                if(ts < it->first)
                    return it->first-ts;                                        // ts falls before the first segment
                // otherwise:
                while(it != lines[f].cend())
                {
                    auto tsi = it->first;
                    auto tfi = it->second;
                    if(ts > tsi && ts < tfi) return -(tfi - ts);                // ts falls withing segment
                    ++it;
                    if(ts > tfi)                                                // in free space behind segment
                    {   
                        if(it == lines[f].cend())                              // ts falls behind last segment
                        {
                           it = lines[f].cbegin();
                           return (it->first + conf.beta - ts);
                        }
                        else if(ts < it->first)                                // ts falls in free space
                            return (it->first - ts);
                    }
                }
                LOG(WARNING) << "Maximum length could not be computed. ts == ts1 in list. This is unlikely.";
                return 0;
            }


            int insertSegment(const RealT ts, const RealT tf, const unsigned char f)
            {
                if(maxl(ts, f) < (tf-ts) || tf < ts) return -1;                   // no valid sgment
                auto res = lines[f].insert(std::make_pair(ts, tf));
                if(res.second)
                    return std::distance(lines[f].begin(), res.first);
                return -1;
            }

            SConfig deleteSegment(const unsigned int i, const unsigned char f)
            {
                if(i >= lines[f].size()) return std::make_pair(-1.,-1.);
                auto it = lines[f].begin();
                std::advance (it,i);
                auto c = *it;
                lines[f].erase(it);
                return c;
            }


            SConfig deleteSegmentT(const RealT ts, const unsigned char f)
            { 
                auto it = lines[f].begin();
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
                if(ts <= it->second - conf.beta)
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
                SConfig old = deleteSegmentT(ts, f);
                if(old.first > 0)
                {
                    ts = ts > old.first ? ts : ts + conf.beta;
                    if(!insertSegment(old.first, ts, f))                         // roll back on failure
                    {
                        LOG(WARNING) << "Error during insertion of first segment in anti segment insertion. ts: " << old.first << " - tf: " << ts;
                        insertSegment(old.first, old.second, f);
                        return false;
                    }
                    tf = tf > conf.beta ? tf - conf.beta : tf;
                    RealT tf1 = old.second;
                    if(tf1 - tf > conf.beta)
                        tf1 = old.second - conf.beta;
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
                {
                    insertSegment(0.,conf.beta-std::numeric_limits<RealT>::epsilon(), f);
                    return true;
                }
                // more than one segment present, continue normally
                SConfig old1 = deleteSegment( (i > 0 ? i-1 : lines[f].size()-1), f) ;
                if(old1.first < 0)
                {
                    insertSegment(old2.first, old2.second, f);
                    LOG(WARNING) << "Could not remove first old segment for anti segment removal. ts: " << old2.first << " i: " << i;
                    return false;
                }
                LOG(INFO) << " :: " << old1.first << " ---- " << old2.second;
                RealT tf = old1.first < old2.second ? old2.second : old2.second + conf.beta;
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
                    for(auto it = lines[fi].begin(); it != lines[fi].end(); it++)
                    {
                        tmp[fi] += " -- " + std::to_string(it->first) + "=="  + std::to_string(it->second);
                    }
                    res += tmp[fi] + "\n";
                }
                return res;
            }
            
        private:
            const Config &conf;
            typedef std::set<SConfig, compare> SConfigL;                                  // Time orderer list (set)
            SConfigL lines[NFLAVORS];
    }; 
}
#endif
