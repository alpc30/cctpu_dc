#ifndef UPDATEQUEUE_H_
#define UPDATEQUEUE_H_

//#pragma once
#include<deque>

namespace TRANSPLAN
{
  
  template <class T> 
    class UpdateQueue
    {
    public:
      std::deque<T> list;

      UpdateQueue()
      {
	
      }
      
      UpdateQueue(const UpdateQueue& queue)
      :list(queue.list)
      {
	
      }
    
      UpdateQueue& operator=(const UpdateQueue& queue)
      {
	this->list = queue.list;
	
	return *this;
      }


      bool empty()
      {
	return list.empty();
      }

      int size()
      {
	return list.size();
      }

      void clear()
      {
	list.clear();
      }

      void push_back(T element);

      T pop_front();

    };
    
 
  template <class T>
    void UpdateQueue<T>::push_back(T element)
    {
      bool exists = false;
      
      
      for (int i = 0; i < list.size(); ++i )
      {
	if ( list[i] == element)
	{
	  exists = true;
	  break;
	}
      }
      
      if ( !exists )
      {
	list.push_back(element);
      }
    }

  template <class T>
    T UpdateQueue<T>::pop_front()
    {
      T element = list[0];
      list.pop_front();
      return element;
    }

}
#endif
