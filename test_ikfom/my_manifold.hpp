
#include "../include/IKFoM_toolkit/mtk/src/mtkmath.hpp"
#include "../include/IKFoM_toolkit/mtk/types/vect.hpp"
#include "../include/IKFoM_toolkit/mtk/types/SOn.hpp"
#include "../include/IKFoM_toolkit/mtk/build_manifold.hpp"

/* this struct defines the compund manifold composed by:
    1. vect3 force robot 1
    2. vect3 torque robot 1
    3. vect3 force robot 2
    4. vect3 torque robot 2
    5. vect3 position robot 1
    6. SO3 orientation robot1
    7. vect3 position robot 2
    8. SO3 orientation robot2
    9. std::vector<vect3,num_pose_measure> position object
    10. std::vector<SO3,num_pose_measure> orientation object
*/

typedef MTK::vect<3, double> vect3;
typedef MTK::SO3<double> SO3;

// define base manifolds
MTK_BUILD_MANIFOLD(pose_measure, ((vect3, position))((SO3, quaternion)));
MTK_BUILD_MANIFOLD(wrench_measure, ((vect3, force))((vect3, tau)));

#undef vect3
#undef SO3

namespace MTK
{
template <class _scalar = double, int num_pose_measure = 1>
struct MyManifold
{
  enum
  {
    DOF = 12 + 12 + 6 * num_pose_measure,
    DIM = 12 + 12 + 6 * num_pose_measure,
    TYP = 5
  };
  typedef _scalar scalar;
  wrench_measure e1he1_e1;
  wrench_measure e2he2_e2;
  pose_measure b1_pose_e1;
  pose_measure b2_pose_e2;
  std::vector<pose_measure> object_pose;

  MyManifold(const wrench_measure& e1he1_e1 = wrench_measure(), const wrench_measure& e2he2_e2 = wrench_measure(),
             const pose_measure& b1_pose_e1 = pose_measure(), const pose_measure& b2_pose_e2 = pose_measure(),
             const std::vector<pose_measure>& object_pose = std::vector<pose_measure>())
    : e1he1_e1(e1he1_e1), e2he2_e2(e2he2_e2), b1_pose_e1(b1_pose_e1), b2_pose_e2(b2_pose_e2)
  {
    this->object_pose.resize(num_pose_measure);
    if (object_pose.size() != num_pose_measure)
    {
      std::cerr << "Error: object_pose size does not match num_pose_measure." << std::endl;
      std::exit(1);
    }
    if (object_pose.empty())
    {
      for (int i = 0; i < num_pose_measure; ++i)
      {
        this->object_pose[i] = pose_measure();
      }
    }
    else if (object_pose.size() != num_pose_measure)
    {
      std::cerr << "Error: object_pose size does not match num_pose_measure." << std::endl;
      std::exit(1);
    }
    for (int i = 0; i < num_pose_measure; ++i)
    {
      this->object_pose[i] = object_pose[i];
    }
  }

  MyManifold()
  {
    this->e1he1_e1 = wrench_measure();
    this->e2he2_e2 = wrench_measure();
    this->b1_pose_e1 = pose_measure();
    this->b2_pose_e2 = pose_measure();
    this->object_pose.resize(num_pose_measure);
    for (int i = 0; i < num_pose_measure; ++i)
    {
      this->object_pose[i] = pose_measure();
    }
  }

  // Manifold requirements

  void boxplus(MTK::vectview<const scalar, DOF> vec, scalar scale = 1)
  {
    e1he1_e1.boxplus(vec.template segment<6>(0), scale);
    e2he2_e2.boxplus(vec.template segment<6>(6), scale);
    b1_pose_e1.boxplus(vec.template segment<6>(12), scale);
    b2_pose_e2.boxplus(vec.template segment<6>(18), scale);
    for (int i = 0; i < num_pose_measure; ++i)
    {
      object_pose[i].boxplus(vec.template segment<6>(24 + 6 * i), scale);
    }
  }
  void boxminus(MTK::vectview<scalar, DOF> res, const MyManifold<scalar, num_pose_measure>& other) const
  {
    e1he1_e1.boxminus(res.template segment<6>(0), other.e1he1_e1);
    e2he2_e2.boxminus(res.template segment<6>(6), other.e2he2_e2);
    b1_pose_e1.boxminus(res.template segment<6>(12), other.b1_pose_e1);
    b2_pose_e2.boxminus(res.template segment<6>(18), other.b2_pose_e2);
    for (int i = 0; i < num_pose_measure; ++i)
    {
      object_pose[i].boxminus(res.template segment<6>(24 + 6 * i), other.object_pose[i]);
    }
  }

  void oplus(MTK::vectview<const scalar, DOF> vec, scalar scale = 1)
  {
    e1he1_e1.oplus(vec.template segment<6>(0), scale);
    e2he2_e2.oplus(vec.template segment<6>(6), scale);
    b1_pose_e1.oplus(vec.template segment<6>(12), scale);
    b2_pose_e2.oplus(vec.template segment<6>(18), scale);
    for (int i = 0; i < num_pose_measure; ++i)
    {
      object_pose[i].oplus(vec.template segment<6>(24 + 6 * i), scale);
    }
  }
};

}  // namespace MTK
