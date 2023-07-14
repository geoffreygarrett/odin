#ifndef ODIN_TRAVERSING_HPP
#define ODIN_TRAVERSING_HPP

#include "include/odin/tree/base.hpp"
#include <functional>
#include <mutex>
#include <odin/tree/base.hpp>

namespace odin {

    template<class Derived, class T>
    class BaseTraverser {
    public:
        using tree_type = Tree<T>;
        using node_type = typename tree_type::raw_node_type;

        explicit BaseTraverser(tree_type &tree) : m_tree(tree) {}

    protected:
        tree_type &m_tree;
    };

    template<class T>
    class SerialTraverser : public BaseTraverser<SerialTraverser<T>, T> {
    public:
        using BaseTraverser<SerialTraverser<T>, T>::BaseTraverser;
        using tree_type = Tree<T>;
        using node_type = typename tree_type::raw_node_type;

        void pre_order(std::function<void(const T &)> visit) const {
            // implement pre-order traversal using m_tree
        }

        void pre_order(std::function<void(T &)> visit) {
            // implement pre-order traversal using m_tree
        }

        void post_order(std::function<void(const T &)> visit) const {
            // implement post-order traversal using m_tree
        }

        void post_order(std::function<void(T &)> visit) {
            // implement post-order traversal using m_tree
        }
    };

    template<class T>
    class ParallelTraverser : public BaseTraverser<ParallelTraverser<T>, T> {
    public:
        using BaseTraverser<ParallelTraverser<T>, T>::BaseTraverser;
        using tree_type = Tree<T>;
        using node_type = typename tree_type::raw_node_type;

        void pre_order(std::function<void(const T &)> visit) const {
            // implement parallel pre-order traversal using m_tree
        }

        void pre_order(std::function<void(T &)> visit) {
            // implement parallel pre-order traversal using m_tree
        }

        void post_order(std::function<void(const T &)> visit) const {
            // implement parallel post-order traversal using m_tree
        }

        void post_order(std::function<void(T &)> visit) {
            // implement parallel post-order traversal using m_tree
        }
    };

    template<class T>
    class MutexTraverser : public BaseTraverser<MutexTraverser<T>, T> {
    public:
        using BaseTraverser<MutexTraverser<T>, T>::BaseTraverser;
        using tree_type = Tree<T>;
        using node_type = typename tree_type::raw_node_type;

        void pre_order(std::function<void(T &)> visit) {
            // implement mutex-protected pre-order traversal using m_tree
        }

        void post_order(std::function<void(T &)> visit) {
            // implement mutex-protected post-order traversal using m_tree
        }

    private:
        std::mutex m_mutex;
    };

    template<class T>
    class UnsafeTraverser : public BaseTraverser<UnsafeTraverser<T>, T> {
    public:
        using BaseTraverser<UnsafeTraverser<T>, T>::BaseTraverser;
        using tree_type = Tree<T>;
        using node_type = typename tree_type::raw_node_type;

        void pre_order(std::function<void(T *)> visit) {
            // implement potentially unsafe pre-order traversal using m_tree
        }

        void post_order(std::function<void(T *)> visit) {
            // implement potentially unsafe post-order traversal using m_tree
        }
    };

}// namespace odin

#endif//ODIN_TRAVERSING_HPP

//
//// Overloaded const-correct serial pre-order traversal
//void pre_order_traversal(std::function<void(const T &)> visit) const {
//    std::function<void(const raw_node_type *)> worker
//            = [&visit, &worker](const raw_node_type *node) {
//                  if (node) {
//                      visit(node->get_data());
//                      for (auto &child: node->get_children()) { worker(child.get()); }
//                  }
//              };
//    worker(m_root.get());
//}
//
//// Overloaded potentially unsafe serial pre-order traversal
//void pre_order_traversal(std::function<void(T *)> visit) {
//    std::function<void(raw_node_type *)> worker = [&visit, &worker](raw_node_type *node) {
//        if (node) {
//            visit(node->get_data_ptr());
//            for (auto &child: node->get_children()) { worker(child.get()); }
//        }
//    };
//    worker(m_root.get());
//}
//
//
//// Parallel Pre-order traversal
//void pre_order_traversal_mt(std::function<void(T *)> visit) {
//    tbb::task_group group;
//
//    std::function<void(raw_node_type *)> worker
//            = [&visit, &group, &worker](raw_node_type *node) {
//                  if (node) {
//                      visit(node->get_data_ptr());
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                  }
//              };
//    worker(m_root.get());
//    group.wait();
//}
//
//// Thread-safe Pre-order traversal (Multi-threaded mutex)
//void pre_order_traversal_mtm(std::function<void(T &)> visit) {
//    tbb::task_group group;
//    std::mutex      mutex;
//
//    std::function<void(raw_node_type *)> worker
//            = [&visit, &group, &worker, &mutex](raw_node_type *node) {
//                  if (node) {
//                      std::lock_guard<std::mutex> lock(mutex);
//                      visit(node->get_data());
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                  }
//              };
//
//    worker(m_root.get());
//
//    group.wait();
//}
//
//// Serial Post-order traversal
//void post_order_traversal(std::function<void(T &)> visit) {
//    std::function<void(raw_node_type *)> worker = [&visit, &worker](raw_node_type *node) {
//        if (node) {
//            for (auto &child: node->get_children()) { worker(child.get()); }
//            visit(node->get_data());
//        }
//    };
//
//    worker(m_root.get());
//}
//
//// Parallel Pre-order traversal (Thread-Safe)
//void pre_order_traversal_mt(const std::function<void(const T &)> &visit) const {
//    tbb::task_group                            group;
//    std::function<void(const raw_node_type *)> worker
//            = [&visit, &group, &worker](const raw_node_type *node) {
//                  if (node) {
//                      visit(node->get_data());
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                  }
//              };
//    worker(m_root.get());
//    group.wait();
//}
//
//// Overloaded potentially unsafe Parallel Pre-order traversal
//void pre_order_traversal_mt(const std::function<void(T *)> &visit) {
//    tbb::task_group                      group;
//    std::function<void(raw_node_type *)> worker
//            = [&visit, &group, &worker](raw_node_type *node) {
//                  if (node) {
//                      visit(node->get_data_ptr());
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                  }
//              };
//    worker(m_root.get());
//    group.wait();
//}
//
//// Parallel Post-order traversal (Thread-Safe)
//void post_order_traversal_mt(const std::function<void(const T &)> &visit) const {
//    tbb::task_group                            group;
//    std::function<void(const raw_node_type *)> worker
//            = [&visit, &group, &worker](const raw_node_type *node) {
//                  if (node) {
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                      visit(node->get_data());
//                  }
//              };
//    worker(m_root.get());
//    group.wait();
//}
//
//// Overloaded potentially unsafe Parallel Post-order traversal
//void post_order_traversal_mt(const std::function<void(T *)> &visit) {
//    tbb::task_group                      group;
//    std::function<void(raw_node_type *)> worker
//            = [&visit, &group, &worker](raw_node_type *node) {
//                  if (node) {
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                      visit(node->get_data_ptr());
//                  }
//              };
//    worker(m_root.get());
//    group.wait();
//}
//
//// Thread-safe Post-order traversal
//void post_order_traversal_mtm(std::function<void(T &)> visit) {
//    tbb::task_group group;
//    std::mutex      mutex;
//
//    std::function<void(raw_node_type *)> worker
//            = [&visit, &group, &worker, &mutex](raw_node_type *node) {
//                  if (node) {
//                      for (auto &child: node->get_children()) {
//                          group.run([&]() { worker(child.get()); });
//                      }
//                      std::lock_guard<std::mutex> lock(mutex);
//                      visit(node->get_data());
//                  }
//              };
//
//    worker(m_root.get());
//
//    group.wait();
//}