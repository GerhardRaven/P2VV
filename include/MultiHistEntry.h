// $Id: $
#ifndef MULTIHISTENTRY_H 
#define MULTIHISTENTRY_H 1

// Include files
#include <RooCategoryProxy.h>

/** @class MultiHistEntry MultiHistEntry.h
 *  
 *
 *  @author Roel Aaij
 *  @date   2012-08-21
 */
template <class EFFICIENCY, class PARENT> 
class MultiHistEntry {
public:

   MultiHistEntry()
   : m_rawEff(0), m_rawRel(0), m_efficiency(0), m_relative(0),
        m_index(0)
   {
   }

   MultiHistEntry(const std::map<RooAbsCategory*, std::string>& categories,
                  EFFICIENCY* efficiency, RooAbsReal* relative)
   : m_rawCats(categories), m_rawEff(efficiency), m_rawRel(relative),
        m_efficiency(0), m_relative(0), m_index(0)
   {
   }

   MultiHistEntry(const MultiHistEntry& other, PARENT* parent)
      : m_rawCats(other.m_rawCats), m_rawEff(other.m_rawEff), m_rawRel(other.m_rawRel),
        m_index(other.m_index)
   {
      if (!other.m_efficiency) {
         m_efficiency = 0;
         m_relative = 0;
         return;
      }
      m_efficiency = new RooRealProxy(other.m_efficiency->GetName(), parent, *other.m_efficiency);
      m_relative = new RooRealProxy(other.m_relative->GetName(), parent, *other.m_relative);
      
      for (std::map<RooCategoryProxy*, std::string>::const_iterator it = other.m_categories.begin(),
              end = other.m_categories.end(); it != end; ++it) {
         m_categories.insert(make_pair(new RooCategoryProxy(it->first->GetName(), parent, *(it->first)),
                                       it->second));
      }
   }

   virtual ~MultiHistEntry()
   {
      for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
              end = m_categories.end(); it != end; ++it) {
         if (it->first) delete it->first;
      }
      m_categories.clear();
      if (m_efficiency) delete m_efficiency;
      if (m_relative) delete m_relative;
   }


   const EFFICIENCY* efficiency() const {
      return const_cast<MultiHistEntry*>(this)->efficiency();
      // }
   }

   EFFICIENCY* efficiency() {
      if (m_efficiency) {
         return dynamic_cast<EFFICIENCY*>(m_efficiency->absArg());
      } else {
         return m_rawEff;
      }
   }

   void setEfficiency(EFFICIENCY* eff) {
      if (m_efficiency) {
         m_efficiency->setArg(*eff);
      } else {
         m_rawEff = eff;
      }
   }

   void setRelative(RooAbsReal* rel) {
      if (m_relative) {
         m_relative->setArg(*rel);
      } else {
         m_rawRel = rel;
      }
   }
   
   RooAbsReal* relative() {
      return m_relative ? dynamic_cast<RooAbsReal*>(m_relative->absArg()) : m_rawRel;
   }

   const RooAbsReal* relative() const{
      return const_cast<MultiHistEntry*>(this)->relative();      
   }

   void setParent(PARENT* parent)
   {
      assert(m_efficiency == 0);
      assert(m_relative == 0);
      
      std::string name;
      for (std::map<RooAbsCategory*, std::string>::const_iterator it = m_rawCats.begin(),
              end = m_rawCats.end(); it != end; ++it) {
         name = it->first->GetName(); name += "_proxy";
         RooCategoryProxy* proxy = new RooCategoryProxy(name.c_str(), name.c_str(), parent,
                                                        *(it->first));
         m_categories.insert(make_pair(proxy, it->second));
      }
      m_rawCats.clear();
      
      name = m_rawEff->GetName(); name += "_proxy";
      m_efficiency = new RooRealProxy(name.c_str(), name.c_str(), parent, *m_rawEff);
      name = m_rawEff->GetName(); name += "_proxy";
      m_relative = new RooRealProxy(name.c_str(), name.c_str(), parent, *m_rawRel);
      
      // m_rawEff = 0;
      // m_rawRel = 0;
   }

   RooArgSet categories() const
   {
      RooArgSet r;
      if (!m_rawCats.empty()) {
         for (std::map<RooAbsCategory*, std::string>::const_iterator it = m_rawCats.begin(),
                 end = m_rawCats.end(); it != end; ++it) {
            if (!it->first) continue;
            r.add(*(it->first));
         }
      } else {
         for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
                 end = m_categories.end(); it != end; ++it) {
            if (!it->first) continue;
            r.add(it->first->arg());
         }
      }
      return r;
   }

   bool thisEntry() const
   {
      bool r = true;
      for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
              end = m_categories.end(); it != end; ++it) {
         const RooAbsCategory* cat = static_cast<const RooAbsCategory*>(it->first->absArg());
         if (strcmp(cat->getLabel(), it->second.c_str()) != 0) {
            r = false;
            break;
         }
      }
      return r;
   }

   void setIndex(const Int_t index) {
      m_index = index;
   }
   
   Int_t index() const {
      return m_index;
   }

   void select()
   {
      if (!m_rawCats.empty()) {
         for (std::map<RooAbsCategory*, std::string>::const_iterator it = m_rawCats.begin(),
                 end = m_rawCats.end(); it != end; ++it) {
            if (!it->first) continue;
            RooAbsCategoryLValue* lval = dynamic_cast<RooAbsCategoryLValue*>(it->first);
            lval->setLabel(it->second.c_str());
         }
      } else {
         for (std::map<RooCategoryProxy*, std::string>::const_iterator it = m_categories.begin(),
                 end = m_categories.end(); it != end; ++it) {
            RooAbsCategoryLValue* lval = dynamic_cast<RooAbsCategoryLValue*>(it->first->absArg());
            lval->setLabel(it->second.c_str());
         }
      }
   }

private:


   std::map<RooAbsCategory*, std::string> m_rawCats;
   EFFICIENCY* m_rawEff; //!
   RooAbsReal* m_rawRel; //!

   std::map<RooCategoryProxy*, std::string> m_categories;
   RooRealProxy* m_efficiency;
   RooRealProxy* m_relative;
   Int_t m_index;

};
#endif // MULTIHISTENTRY_H
