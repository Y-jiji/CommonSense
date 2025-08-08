/* This is wrapper for a priority queue with updatable priorities.
There are three possible implementations.*/

#pragma once

#include <exception>
#include <functional>
#include <queue>
#include <set>
#include <unordered_map>

#include "updatable_priority_queue/updatable_priority_queue.h"

namespace ONIAK
{

  enum class UpdatePQBackend
  {
    PriorityQueue,
    Set,
    UpdatablePriorityQueue,
  };

  // customers with higher priorities are served first
  template <typename Key, typename Val>
  class PriorityUpdatePQ
  {
  public:
    using VKPair = std::pair<Val, Key>;

    PriorityUpdatePQ(std::function<Val(Val)> func = std::identity::operator()<Val>) : val_to_priority_(func) {}

    void push(const Key &key, const Val &val)
    {
      kvmap_[key] = val;
      pq_.emplace(val_to_priority_(val), key);
    }

    void update(const Key &key, const Val &val)
    {
      push(key, val);
    }

    VKPair top() const
    {
      if (!empty())
      {
        auto [val, key] = pq_.top();
        return {kvmap_.at(key), key};
      }
      else
        throw std::runtime_error("Priority queue is empty!");
    }

    void pop()
    {
      if (!empty())
      {
        //   auto [val, key] = pq_.top();
        pq_.pop();
        //   kvmap_.erase(key); keeps consistency of values
      }
      else
        throw std::runtime_error("Priority queue is empty!");
    }

    bool empty() const
    {
      while (!pq_.empty())
      {
        VKPair result = pq_.top();
        auto [val, key] = result;
        if (kvmap_.contains(key) && val == val_to_priority_(kvmap_.at(key)))
          return false;
        else
          pq_.pop();
      }
      return true;
    }

    // This only sets kvmap. Does not push the priority queue.
    void set(const Key &key, const Val &val)
    {
      kvmap_[key] = val;
    }

    const Val &operator[](const Key &key)
    {
      return kvmap_[key];
    }

    void clear()
    {
      pq_ = {};
      kvmap_.clear();
    }

    int size() const
    {
      return pq_.size();
    }

    const auto kvmap() const
    {
      return kvmap_;
    }

  private:
    mutable std::priority_queue<VKPair> pq_;
    std::unordered_map<Key, Val> kvmap_;
    // function that maps val to priority. Default is identity.
    std::function<Val(Val)> val_to_priority_;
  };

  template <typename Key, typename Val>
  class SetPQ
  {
  public:
    using VKPair = std::pair<Val, Key>;
    SetPQ(std::function<Val(Val)> func = std::identity::operator()<Val>) : val_to_priority_(func) {}

    void push(const Key &key, const Val &val)
    {
      kvmap_[key] = val;
      set_.emplace(val_to_priority_(val), key);
    }

    void update(const Key &key, const Val &val)
    {
      set_.erase({val_to_priority_(kvmap_[key]), key});
      kvmap_[key] = val;
      set_.emplace(val_to_priority_(val), key);
    }

    VKPair top() const
    {
      if (!empty())
      {
        auto [priority, key] = *set_.begin();
        return {kvmap_.at(key), key};
      }
      else
        throw std::runtime_error("Priority queue is empty!");
    }

    void pop()
    {
      if (!empty())
      {
        auto [priority, key] = *set_.begin();
        set_.erase({priority, key});
        kvmap_.erase(key);
      }
      else
        throw std::runtime_error("Priority queue is empty!");
    }

    bool empty() const
    {
      return set_.empty();
    }

    // Erases original vkpair from set.
    void set(const Key &key, const Val &val)
    {
      set_.erase({val_to_priority_(kvmap_[key]), key});
      kvmap_[key] = val;
    }

    // This only sets kvmap. Does not push the priority queue.
    const Val &operator[](const Key &key)
    {
      return kvmap_[key];
    }

    void clear()
    {
      set_.clear();
      kvmap_.clear();
    }

    int size() const
    {
      return set_.size();
    }

    const auto kvmap() const
    {
      return kvmap_;
    }

    bool contains(const Key &key)
    {
      return set_.contains({val_to_priority_(kvmap_[key]), key});
    }

  private:
    std::set<VKPair, std::greater<VKPair>> set_;
    std::unordered_map<Key, Val> kvmap_;
    // function that maps val to priority. Default is identity.
    std::function<Val(Val)> val_to_priority_;
  };

  // customers with higher priorities are served first

  // known issue: if val_to_priority_ is not a bijective function,
  // say with f(a) = f(b).
  // If an element was inserted with priority a, and then its priority
  // is updated to b (no insertion), then an element with priority b will 
  // still be popped (whereas it should not be there).

  template <typename Key, typename Val>
  class BetterUpdatePQ
  {
  public:
    using VKPair = std::pair<Val, Key>;

    BetterUpdatePQ(std::function<Val(Val)> func = std::identity::operator()<Val>) : val_to_priority_(func) {}

    void push(const Key &key, const Val &val)
    {
      kvmap_[key] = val;
      pq_.push(key, val_to_priority_(val));
    }

    void update(const Key &key, const Val &val)
    {
      kvmap_[key] = val;
      pq_.set(key, val_to_priority_(val));
    }

    VKPair top() const
    {
      if (!empty())
      {
        auto [val, key] = pq_.top();
        return {kvmap_.at(key), key};
      }
      else
        throw std::runtime_error("Priority queue is empty!");
    }

    void pop()
    {
      if (!empty())
      {
        auto [val, key] = pq_.top();
        pq_.pop();
        kvmap_.erase(key);
      }
      else
        throw std::runtime_error("Priority queue is empty!");
    }

    bool empty() const
    {
      return pq_.empty(); // || pq_.top().priority == -1;
    }

    bool contains(const Key &key)
    {
      return pq_.contains(key);
    }

    // This only sets kvmap. Removes previous value from queue
    void set(const Key &key, const Val &val)
    {
      kvmap_[key] = val;
      pq_.erase(key);
    }

    const Val &operator[](const Key &key)
    {
      return kvmap_[key];
    }

    void clear()
    {
      pq_ = {};
      kvmap_ = {};
    }

    int size() const
    {
      return pq_.size();
    }

    const auto kvmap() const
    {
      return kvmap_;
    }

  private:
    better_priority_queue::updatable_priority_queue<Key, Val> pq_;
    std::unordered_map<Key, Val> kvmap_;
    // function that maps val to priority. Default is identity.
    std::function<Val(Val)> val_to_priority_;
  };

  template <typename Key, typename Val, UpdatePQBackend backend>
  class UpdatePQAdapter
  {
  public:
    using type = void;
    UpdatePQAdapter()
    {
      static_assert(backend == UpdatePQBackend::PriorityQueue, "Backend not supported!");
      // If backend is priority queue, the next template will be instantiated.
      // So here static assert will always trigger.
    }
  };

  // In this implementation, elements are not deleted, only marked as invalid.
  // A key-to-value map is used to keep track of the latest version.
  template <typename Key, typename Val>
  class UpdatePQAdapter<Key, Val, UpdatePQBackend::PriorityQueue>
  {
  public:
    using type = PriorityUpdatePQ<Key, Val>;
  };

  // implementation from github.com/Ten0/updatable_priority_queue
  template <typename Key, typename Val>
  class UpdatePQAdapter<Key, Val, UpdatePQBackend::UpdatablePriorityQueue>
  {
  public:
    using type = BetterUpdatePQ<Key, Val>;
  };

  // implementation from set, priority enforced at insertion
  template <typename Key, typename Val>
  class UpdatePQAdapter<Key, Val, UpdatePQBackend::Set>
  {
  public:
    using type = SetPQ<Key, Val>;
  };

} // namespace ONIAK
