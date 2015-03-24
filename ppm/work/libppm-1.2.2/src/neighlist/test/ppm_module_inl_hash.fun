test_suite ppm_module_inl_hash

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

type(ppm_htable),pointer   :: ht
integer            :: info
integer, parameter :: size = 200

  setup
  nullify(ht)
  allocate(ht)
  end setup

  teardown
  if (associated(ht)) then
  deallocate(ht)
  nullify(ht)
  endif
  end teardown


  test basic
    ! create
    call create_htable(ht, size, info)
    assert_true(info .eq. 0)

    ! basic usage
    assert_true(hash_search(ht, 7_8) .eq. htable_null)
    call hash_insert(ht, 12_8, 144, info)
    assert_true(info .eq. 0)
    assert_true(hash_search(ht, 12_8) .eq. 144)

    assert_true(hash_search(ht, 7_8) .eq. htable_null)
    call hash_insert(ht, 7_8, 49, info)
    assert_true(info .eq. 0)
    assert_true(hash_search(ht, 7_8) .eq. 49)

    ! destroy
    call destroy_htable(ht, info)
    assert_true(info .eq. 0)
  end test


end test_suite
